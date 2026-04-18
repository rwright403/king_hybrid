from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, root_scalar
import CoolProp.CoolProp as CP

from src.models._thermo.n2o_thermo_span_wagner_class import (
    SpanWagnerEOS_SingleState,
    SpanWagnerEOS_EquilibriumPhase,
    lightweight_span_wagner_eos_pressure,
)

# If you prefer to keep your existing HT correlations, import them here instead.
# These names match the ones you have already used in your repo.
from src.models._thermo.convection_heat_transfer import (
    solve_Q_dot_natural_convection_gas,
    solve_Q_dot_natural_convection_liq,
)
from src.models._thermo.conduction_heat_transfer import solve_Q_dot_conduction


R_U = 8.31446261815324
MW_N2O = 44.0128e-3
CW_AL = 896.0


@dataclass
class TankGeometry:
    V_tank: float
    diam_in: float
    diam_out: float

    def __post_init__(self) -> None:
        if self.diam_in >= self.diam_out:
            raise ValueError("diam_in must be smaller than diam_out")
        self.A_cross = 0.25 * np.pi * self.diam_in**2
        self.A_wall_cross = 0.25 * np.pi * (self.diam_out**2 - self.diam_in**2)
        self.height = self.V_tank / self.A_cross
        self.endcap_area = self.A_cross


@dataclass
class Environment:
    T_atm: float
    P_atm: float
    rho_atm: float
    rho_wall: float
    k_wall: float


@dataclass
class ModelParams:
    E_blowing: float = 2.1e4
    tau_cond: float = 0.10
    m_min: float = 1e-9
    rho_min: float = 1e-6
    inner_xtol: float = 1e-10
    step_xtol: float = 1e-10
    step_ftol: float = 1e-10
    max_nfev: int = 100
    gas_T_max: float = 1000.0


@dataclass
class TankState:
    T_liq: float
    T_gas: float
    m_liq: float
    m_gas: float
    T_wall_liq: float
    T_wall_gas: float
    rho_liq: float
    rho_gas: float
    P_tank: float


@dataclass
class TankHistory:
    time: list[float] = field(default_factory=list)
    P_tank: list[float] = field(default_factory=list)
    m_liq: list[float] = field(default_factory=list)
    m_gas: list[float] = field(default_factory=list)
    m_total: list[float] = field(default_factory=list)
    T_liq: list[float] = field(default_factory=list)
    T_gas: list[float] = field(default_factory=list)
    T_sat: list[float] = field(default_factory=list)
    T_wall_liq: list[float] = field(default_factory=list)
    T_wall_gas: list[float] = field(default_factory=list)
    fill_fraction: list[float] = field(default_factory=list)
    m_dot_evap: list[float] = field(default_factory=list)
    m_dot_cond: list[float] = field(default_factory=list)
    m_dot_liq_in: list[float] = field(default_factory=list)
    m_dot_inj_out: list[float] = field(default_factory=list)
    m_dot_vent: list[float] = field(default_factory=list)
    energy_residual: list[float] = field(default_factory=list)


@dataclass
class PublishedReference:
    time: np.ndarray
    P_tank: Optional[np.ndarray] = None
    T_liq: Optional[np.ndarray] = None
    T_gas: Optional[np.ndarray] = None
    T_sat: Optional[np.ndarray] = None
    m_liq: Optional[np.ndarray] = None
    m_gas: Optional[np.ndarray] = None


class ImplicitTankModel:
    """
    Sign convention:
      m_dot_liq_in > 0   entering liquid node
      m_dot_inj_out > 0  liquid leaving tank downstream
      m_dot_vent > 0     gas leaving tank
      m_dot_evap > 0     liquid -> gas
      m_dot_cond > 0     gas -> liquid
    """

    def __init__(
        self,
        name: str,
        geom: TankGeometry,
        env: Environment,
        params: ModelParams,
        initial_state: TankState,
    ) -> None:
        self.name = name
        self.geom = geom
        self.env = env
        self.params = params
        self.state = initial_state
        self.last_diag: Dict[str, Any] = {}

    def _wall_split(self, V_gas: float) -> Dict[str, float]:
        h_gas = np.clip(V_gas / self.geom.A_cross, 0.0, self.geom.height)
        h_liq = np.clip(self.geom.height - h_gas, 0.0, self.geom.height)

        A_liq_wet = np.pi * self.geom.diam_in * h_liq + self.geom.endcap_area
        A_gas_wet = np.pi * self.geom.diam_in * h_gas + self.geom.endcap_area
        A_liq_outer = np.pi * self.geom.diam_out * h_liq + self.geom.endcap_area
        A_gas_outer = np.pi * self.geom.diam_out * h_gas + self.geom.endcap_area

        m_wall_liq = self.env.rho_wall * self.geom.A_wall_cross * h_liq
        m_wall_gas = self.env.rho_wall * self.geom.A_wall_cross * h_gas

        return {
            "h_gas": h_gas,
            "h_liq": h_liq,
            "A_liq_wet": A_liq_wet,
            "A_gas_wet": A_gas_wet,
            "A_liq_outer": A_liq_outer,
            "A_gas_outer": A_gas_outer,
            "m_wall_liq": m_wall_liq,
            "m_wall_gas": m_wall_gas,
        }

    def _thermo_closure(
        self,
        T_liq: float,
        T_gas: float,
        m_liq: float,
        m_gas: float,
        rho_guess: np.ndarray,
    ):
        m_liq = max(m_liq, self.params.m_min)
        m_gas = max(m_gas, self.params.m_min)

        def residual(rhos: np.ndarray) -> np.ndarray:
            rho_liq = max(rhos[0], self.params.rho_min)
            rho_gas = max(rhos[1], self.params.rho_min)
            P_liq = lightweight_span_wagner_eos_pressure(rho_liq, T_liq)
            P_gas = lightweight_span_wagner_eos_pressure(rho_gas, T_gas)
            V_err = (m_liq / rho_liq) + (m_gas / rho_gas) - self.geom.V_tank
            return np.array([P_liq - P_gas, V_err])

        sol = least_squares(
            residual,
            x0=np.maximum(rho_guess, self.params.rho_min),
            bounds=(
                np.array([self.params.rho_min, self.params.rho_min]),
                np.array([1.3e3, 5.0e2]),
            ),
            x_scale=np.maximum(np.abs(rho_guess), 1.0),
            xtol=self.params.inner_xtol,
            ftol=self.params.inner_xtol,
            gtol=self.params.inner_xtol,
            max_nfev=200,
        )
        if not sol.success:
            raise RuntimeError(f"[{self.name}] thermo closure failed: {sol.message}")

        rho_liq = float(sol.x[0])
        rho_gas = float(sol.x[1])
        P_tank = float(lightweight_span_wagner_eos_pressure(rho_gas, T_gas))

        liq = SpanWagnerEOS_SingleState(rho_liq, T_liq)
        gas = SpanWagnerEOS_SingleState(rho_gas, T_gas)
        sat = SpanWagnerEOS_EquilibriumPhase(None, P_tank)

        V_liq = m_liq / rho_liq
        V_gas = m_gas / rho_gas
        return rho_liq, rho_gas, P_tank, V_liq, V_gas, liq, gas, sat

    def _p_sat_from_T(self, T: float) -> float:
        return float(CP.PropsSI("P", "T", T, "Q", 1, "N2O"))

    def _compute_fluxes(
        self,
        T_liq: float,
        T_gas: float,
        T_wall_liq: float,
        T_wall_gas: float,
        rho_liq: float,
        rho_gas: float,
        P_tank: float,
        V_liq: float,
        V_gas: float,
        liq: Any,
        gas: Any,
        sat: Any,
        inputs: Dict[str, float],
        V_liq_old: float,
        dt: float,
    ) -> Dict[str, float]:
        split = self._wall_split(V_gas)

        T_film_gas = 0.5 * (sat.T + T_gas)
        T_film_liq = 0.5 * (sat.T + T_liq)

        Q_dot_sat_to_gas = solve_Q_dot_natural_convection_gas(
            rho_gas, sat.T, T_gas, T_film_gas, P_tank,
            0.15, 0.333, self.geom.diam_in, self.geom.A_cross, "N2O"
        )
        Q_dot_liq_to_sat = self.params.E_blowing * solve_Q_dot_natural_convection_liq(
            rho_liq, T_liq, sat.T, T_film_liq, P_tank,
            0.15, 0.333, self.geom.diam_in, self.geom.A_cross, "N2O"
        )

        delta_h_evap = (sat.h_sat_gas - sat.h_sat_liq) + (sat.h_sat_liq - liq.h)
        if (Q_dot_liq_to_sat - Q_dot_sat_to_gas) > 0.0 and delta_h_evap > 0.0:
            m_dot_evap = (Q_dot_liq_to_sat - Q_dot_sat_to_gas) / delta_h_evap
        else:
            m_dot_evap = 0.0

        try:
            p_sat_gas = self._p_sat_from_T(T_gas)
        except Exception:
            p_sat_gas = np.inf

        if P_tank > p_sat_gas and np.isfinite(p_sat_gas):
            m_dot_cond = ((P_tank - p_sat_gas) * V_gas * MW_N2O) / (
                R_U * max(T_gas, 1.0) * self.params.tau_cond
            )
        else:
            m_dot_cond = 0.0

        Q_dot_liq_wall_to_liq = solve_Q_dot_natural_convection_liq(
            rho_liq, T_wall_liq, T_liq, T_liq, P_tank,
            0.021, 0.4, max(split["h_liq"], 1e-6), split["A_liq_wet"], "N2O"
        )
        Q_dot_gas_wall_to_gas = solve_Q_dot_natural_convection_gas(
            rho_gas, T_wall_gas, T_gas, T_gas, P_tank,
            0.021, 0.4, max(split["h_gas"], 1e-6), split["A_gas_wet"], "N2O"
        )
        Q_dot_atm_to_liq_wall = solve_Q_dot_natural_convection_gas(
            self.env.rho_atm, self.env.T_atm, T_wall_liq, self.env.T_atm, self.env.P_atm,
            0.59, 0.25, self.geom.height, split["A_liq_outer"], "Air"
        )
        Q_dot_atm_to_gas_wall = solve_Q_dot_natural_convection_gas(
            self.env.rho_atm, self.env.T_atm, T_wall_gas, self.env.T_atm, self.env.P_atm,
            0.59, 0.25, self.geom.height, split["A_gas_outer"], "Air"
        )

        # IMPORTANT: your own solve_Q_dot_conduction signature may differ.
        try:
            Q_dot_liq_wall_to_gas_wall = solve_Q_dot_conduction(
                T_wall_liq - T_wall_gas, self.geom.height, self.env.k_wall,
                self.geom.diam_in, self.geom.diam_out
            )
        except TypeError:
            Q_dot_liq_wall_to_gas_wall = solve_Q_dot_conduction(
                T_wall_liq - T_wall_gas, self.geom.height, self.env.k_wall, self.geom.A_wall_cross
            )

        m_dot_liq_in = max(float(inputs.get("m_dot_liq_in", 0.0)), 0.0)
        h_liq_in = float(inputs.get("h_liq_in", liq.h))
        m_dot_inj_out = max(float(inputs.get("m_dot_inj_out", 0.0)), 0.0)
        m_dot_vent = max(float(inputs.get("m_dot_vent", 0.0)), 0.0)

        V_dot_liq = (V_liq - V_liq_old) / dt
        V_dot_gas = -V_dot_liq

        return {
            **split,
            "m_dot_evap": m_dot_evap,
            "m_dot_cond": m_dot_cond,
            "m_dot_liq_in": m_dot_liq_in,
            "h_liq_in": h_liq_in,
            "m_dot_inj_out": m_dot_inj_out,
            "m_dot_vent": m_dot_vent,
            "Q_dot_sat_to_gas": Q_dot_sat_to_gas,
            "Q_dot_liq_to_sat": Q_dot_liq_to_sat,
            "Q_dot_liq_wall_to_liq": Q_dot_liq_wall_to_liq,
            "Q_dot_gas_wall_to_gas": Q_dot_gas_wall_to_gas,
            "Q_dot_atm_to_liq_wall": Q_dot_atm_to_liq_wall,
            "Q_dot_atm_to_gas_wall": Q_dot_atm_to_gas_wall,
            "Q_dot_liq_wall_to_gas_wall": Q_dot_liq_wall_to_gas_wall,
            "Q_dot_liq": Q_dot_liq_wall_to_liq - Q_dot_liq_to_sat,
            "Q_dot_gas": Q_dot_gas_wall_to_gas + Q_dot_sat_to_gas,
            "V_dot_liq": V_dot_liq,
            "V_dot_gas": V_dot_gas,
        }

    def _step_residual(self, y_new: np.ndarray, dt: float, inputs: Dict[str, float]) -> np.ndarray:
        T_liq_n, T_gas_n, m_liq_n, m_gas_n, T_wliq_n, T_wgas_n = y_new

        if m_liq_n <= self.params.m_min or m_gas_n <= self.params.m_min:
            return np.full(6, 1e20)
        if T_liq_n < 150.0 or T_gas_n < 150.0 or T_gas_n > self.params.gas_T_max:
            return np.full(6, 1e20)

        try:
            rho_liq_n, rho_gas_n, P_n, V_liq_n, V_gas_n, liq_n, gas_n, sat_n = self._thermo_closure(
                T_liq_n, T_gas_n, m_liq_n, m_gas_n,
                np.array([self.state.rho_liq, self.state.rho_gas]),
            )

            # Old state is already known; do not re-solve the thermo closure here.
            V_gas_o = self.state.m_gas / self.state.rho_gas
            liq_o = SpanWagnerEOS_SingleState(self.state.rho_liq, self.state.T_liq)
            gas_o = SpanWagnerEOS_SingleState(self.state.rho_gas, self.state.T_gas)
        except Exception:
            return np.full(6, 1e20)

        split_old = self._wall_split(V_gas_o)
        flux = self._compute_fluxes(
            T_liq_n, T_gas_n, T_wliq_n, T_wgas_n,
            rho_liq_n, rho_gas_n, P_n, V_liq_n, V_gas_n,
            liq_n, gas_n, sat_n, inputs, self.state.m_liq / self.state.rho_liq, dt
        )

        m_dot_liq = flux["m_dot_liq_in"] + flux["m_dot_cond"] - flux["m_dot_evap"] - flux["m_dot_inj_out"]
        m_dot_gas = flux["m_dot_evap"] - flux["m_dot_cond"] - flux["m_dot_vent"]

        U_liq_o = self.state.m_liq * liq_o.u
        U_gas_o = self.state.m_gas * gas_o.u
        U_liq_n = m_liq_n * liq_n.u
        U_gas_n = m_gas_n * gas_n.u

        U_wall_liq_o = split_old["m_wall_liq"] * CW_AL * self.state.T_wall_liq
        U_wall_gas_o = split_old["m_wall_gas"] * CW_AL * self.state.T_wall_gas
        U_wall_liq_n = flux["m_wall_liq"] * CW_AL * T_wliq_n
        U_wall_gas_n = flux["m_wall_gas"] * CW_AL * T_wgas_n

        R_m_liq = m_liq_n - self.state.m_liq - dt * m_dot_liq
        R_m_gas = m_gas_n - self.state.m_gas - dt * m_dot_gas

        R_U_liq = U_liq_n - U_liq_o - dt * (
            flux["m_dot_liq_in"] * flux["h_liq_in"]
            + flux["m_dot_cond"] * sat_n.h_sat_liq
            - flux["m_dot_evap"] * liq_n.h
            - flux["m_dot_inj_out"] * liq_n.h
            - P_n * flux["V_dot_liq"]
            + flux["Q_dot_liq"]
        )

        R_U_gas = U_gas_n - U_gas_o - dt * (
            flux["m_dot_evap"] * sat_n.h_sat_gas
            - flux["m_dot_cond"] * gas_n.h
            - flux["m_dot_vent"] * gas_n.h
            - P_n * flux["V_dot_gas"]
            + flux["Q_dot_gas"]
        )

        R_U_wall_liq = U_wall_liq_n - U_wall_liq_o - dt * (
            flux["Q_dot_atm_to_liq_wall"]
            - flux["Q_dot_liq_wall_to_liq"]
            - flux["Q_dot_liq_wall_to_gas_wall"]
        )
        R_U_wall_gas = U_wall_gas_n - U_wall_gas_o - dt * (
            flux["Q_dot_atm_to_gas_wall"]
            - flux["Q_dot_gas_wall_to_gas"]
            + flux["Q_dot_liq_wall_to_gas_wall"]
        )

        raw_resid = np.array([R_m_liq, R_m_gas, R_U_liq, R_U_gas, R_U_wall_liq, R_U_wall_gas], dtype=float)
        scales = np.array([
            max(self.state.m_liq, 1.0),
            max(self.state.m_gas, 1.0),
            max(abs(U_liq_o), 1e3),
            max(abs(U_gas_o), 1e3),
            max(abs(U_wall_liq_o), 1e3),
            max(abs(U_wall_gas_o), 1e3),
        ], dtype=float)
        scaled_resid = raw_resid / scales

        self.last_diag = {
            "rho_liq": rho_liq_n,
            "rho_gas": rho_gas_n,
            "P_tank": P_n,
            "V_liq": V_liq_n,
            "V_gas": V_gas_n,
            "raw_resid": raw_resid,
            "scaled_resid": scaled_resid,
            **flux,
        }

        return scaled_resid

    def step(self, dt: float, inputs: Dict[str, float]) -> TankState:
        prev = TankState(**self.state.__dict__)

        y0 = np.array([
            self.state.T_liq,
            self.state.T_gas,
            self.state.m_liq,
            self.state.m_gas,
            self.state.T_wall_liq,
            self.state.T_wall_gas,
        ], dtype=float)

        lb = np.array([150.0, 150.0, self.params.m_min, self.params.m_min, 150.0, 150.0])
        ub = np.array([400.0, self.params.gas_T_max, 1e6, 1e6, 500.0, 500.0])

        sol = least_squares(
            lambda y: self._step_residual(y, dt, inputs),
            x0=y0,
            bounds=(lb, ub),
            xtol=self.params.step_xtol,
            ftol=self.params.step_ftol,
            gtol=self.params.step_ftol,
            max_nfev=self.params.max_nfev,
        )
        if not sol.success:
            raise RuntimeError(f"[{self.name}] implicit step failed: {sol.message}")

        y = sol.x

        # Re-evaluate the accepted point through the residual function so last_diag is
        # populated for the final solution, but do not call _thermo_closure directly here.
        final_resid = self._step_residual(y, dt, inputs)
        if (not np.all(np.isfinite(final_resid))) or np.max(np.abs(final_resid)) > 1e-3:
            print(f"[{self.name}] least_squares status={sol.status}")
            print(f"[{self.name}] message={sol.message}")
            print(f"[{self.name}] nfev={sol.nfev}")
            print(f"[{self.name}] x={sol.x}")
            print(f"[{self.name}] raw resid={self.last_diag.get('raw_resid')}")
            print(f"[{self.name}] scaled resid={self.last_diag.get('scaled_resid')}")
            print(f"[{self.name}] diag={self.last_diag}")
            raise RuntimeError(f"[{self.name}] accepted step ended at an invalid thermodynamic state")

        rho_liq = float(self.last_diag["rho_liq"])
        rho_gas = float(self.last_diag["rho_gas"])
        P_tank = float(self.last_diag["P_tank"])

        self.state = TankState(
            T_liq=float(y[0]),
            T_gas=float(y[1]),
            m_liq=float(y[2]),
            m_gas=float(y[3]),
            T_wall_liq=float(y[4]),
            T_wall_gas=float(y[5]),
            rho_liq=rho_liq,
            rho_gas=rho_gas,
            P_tank=P_tank,
        )
        self.last_diag["energy_audit"] = self.energy_audit(prev, dt)
        return self.state

    def energy_audit(self, prev_state: TankState, dt: float) -> Dict[str, float]:
        # Do not re-solve thermo closure here. The densities and pressures are already part
        # of the accepted old/new states, so just reconstruct the EOS objects directly.
        V_gas_o = prev_state.m_gas / prev_state.rho_gas
        V_gas_n = self.state.m_gas / self.state.rho_gas
        liq_o = SpanWagnerEOS_SingleState(prev_state.rho_liq, prev_state.T_liq)
        gas_o = SpanWagnerEOS_SingleState(prev_state.rho_gas, prev_state.T_gas)
        liq_n = SpanWagnerEOS_SingleState(self.state.rho_liq, self.state.T_liq)
        gas_n = SpanWagnerEOS_SingleState(self.state.rho_gas, self.state.T_gas)

        split_o = self._wall_split(V_gas_o)
        split_n = self._wall_split(V_gas_n)

        U_old = (
            prev_state.m_liq * liq_o.u
            + prev_state.m_gas * gas_o.u
            + split_o["m_wall_liq"] * CW_AL * prev_state.T_wall_liq
            + split_o["m_wall_gas"] * CW_AL * prev_state.T_wall_gas
        )
        U_new = (
            self.state.m_liq * liq_n.u
            + self.state.m_gas * gas_n.u
            + split_n["m_wall_liq"] * CW_AL * self.state.T_wall_liq
            + split_n["m_wall_gas"] * CW_AL * self.state.T_wall_gas
        )

        rhs_external = dt * (
            self.last_diag["m_dot_liq_in"] * self.last_diag["h_liq_in"]
            - self.last_diag["m_dot_inj_out"] * liq_n.h
            - self.last_diag["m_dot_vent"] * gas_n.h
            + self.last_diag["Q_dot_atm_to_liq_wall"]
            + self.last_diag["Q_dot_atm_to_gas_wall"]
        )

        return {
            "dU_total": U_new - U_old,
            "rhs_external": rhs_external,
            "residual": (U_new - U_old) - rhs_external,
        }

    @property
    def V_liq(self) -> float:
        return self.state.m_liq / self.state.rho_liq

    @property
    def V_gas(self) -> float:
        return self.state.m_gas / self.state.rho_gas

    @property
    def fill_fraction(self) -> float:
        return self.V_liq / self.geom.V_tank

    def current_liq_state(self):
        return SpanWagnerEOS_SingleState(self.state.rho_liq, self.state.T_liq)

    def current_gas_state(self):
        return SpanWagnerEOS_SingleState(self.state.rho_gas, self.state.T_gas)


# --------------------------
# Initialization helpers
# --------------------------

def init_equilibrium_tank_from_pressure(
    name: str,
    m_total: float,
    P_tank: float,
    geom: TankGeometry,
    env: Environment,
    params: ModelParams,
) -> ImplicitTankModel:
    sat = SpanWagnerEOS_EquilibriumPhase(None, P_tank)
    rho_bulk = m_total / geom.V_tank
    x = ((1.0 / rho_bulk) - (1.0 / sat.rho_sat_liq)) / ((1.0 / sat.rho_sat_gas) - (1.0 / sat.rho_sat_liq))
    x = float(np.clip(x, 0.0, 1.0))

    m_gas = x * m_total
    m_liq = (1.0 - x) * m_total

    state = TankState(
        T_liq=sat.T,
        T_gas=sat.T,
        m_liq=m_liq,
        m_gas=m_gas,
        T_wall_liq=env.T_atm,
        T_wall_gas=env.T_atm,
        rho_liq=sat.rho_sat_liq,
        rho_gas=sat.rho_sat_gas,
        P_tank=P_tank,
    )
    return ImplicitTankModel(name, geom, env, params, state)


def init_almost_empty_run_tank(
    name: str,
    P_gas_init: float,
    T_gas_init: float,
    geom: TankGeometry,
    env: Environment,
    params: ModelParams,
    T_liq_seed: Optional[float] = None,
    m_liq_seed: float = 1e-5,
) -> ImplicitTankModel:
    sat = SpanWagnerEOS_EquilibriumPhase(None, max(P_gas_init, 8.9e4))
    T_liq = sat.T if T_liq_seed is None else T_liq_seed
    rho_liq = sat.rho_sat_liq

    V_liq = m_liq_seed / rho_liq
    V_gas = geom.V_tank - V_liq
    if V_gas <= 0.0:
        raise ValueError("liquid seed volume is larger than tank")

    def pressure_error(m_gas: float) -> float:
        rho_gas = m_gas / V_gas
        return lightweight_span_wagner_eos_pressure(rho_gas, T_gas_init) - P_gas_init

    sol = root_scalar(
        pressure_error,
        method="secant",
        x0=0.1,
        x1=0.09,
        xtol=1e-12,
        maxiter=200,
    )
    if not sol.converged:
        raise RuntimeError("could not initialize run tank gas mass")

    m_gas = float(sol.root)
    rho_gas = m_gas / V_gas

    state = TankState(
        T_liq=T_liq,
        T_gas=T_gas_init,
        m_liq=m_liq_seed,
        m_gas=m_gas,
        T_wall_liq=env.T_atm,
        T_wall_gas=env.T_atm,
        rho_liq=rho_liq,
        rho_gas=rho_gas,
        P_tank=P_gas_init,
    )
    return ImplicitTankModel(name, geom, env, params, state)


# --------------------------
# Flow models
# --------------------------

def liquid_orifice_mdot(Cd: float, A: float, rho_up: float, P_up: float, P_down: float) -> float:
    dP = max(P_up - P_down, 0.0)
    if dP <= 0.0 or A <= 0.0 or Cd <= 0.0:
        return 0.0
    return Cd * A * np.sqrt(2.0 * rho_up * dP)


def gas_orifice_mdot(Cd: float, A: float, rho_up: float, P_up: float, P_down: float) -> float:
    dP = max(P_up - P_down, 0.0)
    if dP <= 0.0 or A <= 0.0 or Cd <= 0.0:
        return 0.0
    return Cd * A * np.sqrt(2.0 * rho_up * dP)


# --------------------------
# History / plotting
# --------------------------

def append_history(hist: TankHistory, t: float, tank: ImplicitTankModel) -> None:
    hist.time.append(t)
    hist.P_tank.append(tank.state.P_tank)
    hist.m_liq.append(tank.state.m_liq)
    hist.m_gas.append(tank.state.m_gas)
    hist.m_total.append(tank.state.m_liq + tank.state.m_gas)
    hist.T_liq.append(tank.state.T_liq)
    hist.T_gas.append(tank.state.T_gas)
    try:
        hist.T_sat.append(SpanWagnerEOS_EquilibriumPhase(None, tank.state.P_tank).T)
    except Exception:
        hist.T_sat.append(np.nan)
    hist.T_wall_liq.append(tank.state.T_wall_liq)
    hist.T_wall_gas.append(tank.state.T_wall_gas)
    hist.fill_fraction.append(tank.fill_fraction)
    hist.m_dot_evap.append(float(tank.last_diag.get("m_dot_evap", np.nan)))
    hist.m_dot_cond.append(float(tank.last_diag.get("m_dot_cond", np.nan)))
    hist.m_dot_liq_in.append(float(tank.last_diag.get("m_dot_liq_in", np.nan)))
    hist.m_dot_inj_out.append(float(tank.last_diag.get("m_dot_inj_out", np.nan)))
    hist.m_dot_vent.append(float(tank.last_diag.get("m_dot_vent", np.nan)))
    hist.energy_residual.append(float(tank.last_diag.get("energy_audit", {}).get("residual", np.nan)))


def plot_tank_history(hist: TankHistory, title: str, ref: Optional[PublishedReference] = None) -> None:
    plt.figure(figsize=(15, 5))
    plt.suptitle(title)

    plt.subplot(1, 3, 1)
    plt.plot(hist.time, hist.P_tank, label="model tank")
    if ref is not None and ref.P_tank is not None:
        plt.plot(ref.time, ref.P_tank, label="published tank")
    plt.xlabel("Time (s)")
    plt.ylabel("Pressure (Pa)")
    plt.title("Pressure vs Time")
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(hist.time, hist.m_liq, label="model liquid")
    plt.plot(hist.time, hist.m_gas, label="model gas")
    if ref is not None and ref.m_liq is not None:
        plt.plot(ref.time, ref.m_liq, label="published liquid")
    if ref is not None and ref.m_gas is not None:
        plt.plot(ref.time, ref.m_gas, label="published gas")
    plt.xlabel("Time (s)")
    plt.ylabel("Mass (kg)")
    plt.title("Mass vs Time")
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 3)
    plt.plot(hist.time, hist.T_liq, label="model liquid")
    plt.plot(hist.time, hist.T_gas, label="model gas")
    plt.plot(hist.time, hist.T_sat, label="model sat")
    plt.plot(hist.time, hist.T_wall_liq, label="wall liquid")
    plt.plot(hist.time, hist.T_wall_gas, label="wall gas")
    if ref is not None and ref.T_liq is not None:
        plt.plot(ref.time, ref.T_liq, label="published liquid")
    if ref is not None and ref.T_gas is not None:
        plt.plot(ref.time, ref.T_gas, label="published gas")
    if ref is not None and ref.T_sat is not None:
        plt.plot(ref.time, ref.T_sat, label="published sat")
    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (K)")
    plt.title("Temperature vs Time")
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()


def plot_fill_diagnostics(bot_hist: TankHistory, run_hist: TankHistory) -> None:
    plt.figure(figsize=(15, 8))

    plt.subplot(2, 2, 1)
    plt.plot(run_hist.time, run_hist.fill_fraction)
    plt.xlabel("Time (s)")
    plt.ylabel("Run tank fill fraction")
    plt.grid(True)

    plt.subplot(2, 2, 2)
    plt.plot(bot_hist.time, bot_hist.m_dot_inj_out, label="bottle outflow")
    plt.plot(run_hist.time, run_hist.m_dot_liq_in, label="run tank inflow")
    plt.plot(run_hist.time, run_hist.m_dot_vent, label="run tank vent")
    plt.xlabel("Time (s)")
    plt.ylabel("Mass flow (kg/s)")
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(bot_hist.time, bot_hist.energy_residual, label="bottle")
    plt.plot(run_hist.time, run_hist.energy_residual, label="run tank")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy residual per step (J)")
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 2, 4)
    plt.plot(run_hist.time, run_hist.m_dot_evap, label="evap")
    plt.plot(run_hist.time, run_hist.m_dot_cond, label="cond")
    plt.xlabel("Time (s)")
    plt.ylabel("Interfacial mdot (kg/s)")
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()


# --------------------------
# Coupled bottle -> run tank fill
# --------------------------

def run_coupled_fill(
    bottle: ImplicitTankModel,
    run_tank: ImplicitTankModel,
    dt: float,
    t_end: float,
    target_fill_fraction: float,
    Cd_line: float,
    A_line: float,
    Cd_vent: float,
    A_vent: float,
    verbose_every: int = 100,
):
    t = 0.0
    step_idx = 0

    bottle_hist = TankHistory()
    run_hist = TankHistory()
    append_history(bottle_hist, t, bottle)
    append_history(run_hist, t, run_tank)

    start = time.time()

    while t < t_end and run_tank.fill_fraction < target_fill_fraction:
        bottle_liq = bottle.current_liq_state()
        run_gas = run_tank.current_gas_state()

        m_dot_line = liquid_orifice_mdot(
            Cd_line, A_line,
            bottle.state.rho_liq,
            bottle.state.P_tank,
            run_tank.state.P_tank,
        )

        m_dot_run_vent = gas_orifice_mdot(
            Cd_vent, A_vent,
            run_tank.state.rho_gas,
            run_tank.state.P_tank,
            run_tank.env.P_atm,
        )

        bottle_inputs = {
            "m_dot_liq_in": 0.0,
            "h_liq_in": bottle_liq.h,
            "m_dot_inj_out": m_dot_line,
            "m_dot_vent": 0.0,
        }
        run_inputs = {
            "m_dot_liq_in": m_dot_line,
            "h_liq_in": bottle_liq.h,
            "m_dot_inj_out": 0.0,
            "m_dot_vent": m_dot_run_vent,
        }

        bottle.step(dt, bottle_inputs)
        run_tank.step(dt, run_inputs)

        t += dt
        step_idx += 1

        append_history(bottle_hist, t, bottle)
        append_history(run_hist, t, run_tank)

        if verbose_every > 0 and step_idx % verbose_every == 0:
            print(
                f"t={t:.3f} s | run fill={run_tank.fill_fraction:.4f} | "
                f"P_bot={bottle.state.P_tank/1e5:.2f} bar | "
                f"P_run={run_tank.state.P_tank/1e5:.2f} bar | "
                f"Tg_run={run_tank.state.T_gas:.2f} K"
            )

    runtime = time.time() - start
    print(f"Completed in {runtime:.3f} s wall time, sim t={t:.3f} s")
    return bottle_hist, run_hist


# --------------------------
# Example main
# --------------------------

def main() -> None:
    dt = 1e-4
    t_end = 60.0
    target_fill_fraction = 0.90

    P_atm = 101325.0
    T_atm = 286.5
    rho_atm = 1.225

    env_common = Environment(
        T_atm=T_atm,
        P_atm=P_atm,
        rho_atm=rho_atm,
        rho_wall=2770.0,
        k_wall=237.0,
    )
    params = ModelParams(
        E_blowing=2.1e4,
        tau_cond=0.10,
        gas_T_max=800.0,
        max_nfev=100,
    )

    # Bottle
    bot_geom = TankGeometry(V_tank=0.1, diam_in=0.0254 * 7.0, diam_out=(0.0254 * 7.0) + 2.0 * (0.0254 * 0.25))
    bottle = init_equilibrium_tank_from_pressure(
        name="bottle",
        m_total=29.4835,
        P_tank=45e5,
        geom=bot_geom,
        env=env_common,
        params=params,
    )

    # Run tank
    run_geom = TankGeometry(V_tank=0.0208, diam_in=0.1282, diam_out=0.1413)
    run_tank = init_almost_empty_run_tank(
        name="run_tank",
        P_gas_init=P_atm,
        T_gas_init=T_atm,
        geom=run_geom,
        env=env_common,
        params=params,
    )

    # Line and vent
    Cd_line = 0.50
    A_line = 0.25 * np.pi * (0.0254 * 0.4)**2
    Cd_vent = 0.65
    A_vent = 0.25 * np.pi * (0.0254 * 0.002)**2

    bottle_hist, run_hist = run_coupled_fill(
        bottle=bottle,
        run_tank=run_tank,
        dt=dt,
        t_end=t_end,
        target_fill_fraction=target_fill_fraction,
        Cd_line=Cd_line,
        A_line=A_line,
        Cd_vent=Cd_vent,
        A_vent=A_vent,
        verbose_every=100,
    )

    plot_tank_history(bottle_hist, "Bottle")
    plot_tank_history(run_hist, "Run Tank")
    plot_fill_diagnostics(bottle_hist, run_hist)


if __name__ == "__main__":
    main()
