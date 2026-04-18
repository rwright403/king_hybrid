import os
import numpy as np
from src.models._thermo.n2o_thermo_span_wagner_class import SpanWagnerEOS_EquilibriumPhase

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def build_n2o_rho_sat_liq_table(
    P_min=9.0e4,
    P_max=7.2e6,
    nP=1000,
    out_file=None,
):
    if out_file is None:
        out_file = os.path.join(SCRIPT_DIR, "n2o_rho_sat_liq_table.npz")

    P_grid = np.geomspace(P_min, P_max, nP)   # Pa, log-spaced
    rho_sat_liq = np.empty_like(P_grid)

    for i, P in enumerate(P_grid):
        try:
            SWEOS_eq = SpanWagnerEOS_EquilibriumPhase(None, P)
            rho_sat_liq[i] = SWEOS_eq.rho_sat_liq
        except Exception:
           rho_sat_liq[i] = np.nan

    np.savez(out_file, P_grid=P_grid, rho_sat_liq=rho_sat_liq)
    print("saved to:", out_file)
    return P_grid, rho_sat_liq

def build_n2o_rho_sat_gas_table(
    P_min=9.0e4,
    P_max=7.2e6,
    nP=1000,
    out_file=None,
):
    if out_file is None:
        out_file = os.path.join(SCRIPT_DIR, "n2o_rho_sat_gas_table.npz")

    P_grid = np.geomspace(P_min, P_max, nP)   # Pa, log-spaced
    rho_sat_gas = np.empty_like(P_grid)

    for i, P in enumerate(P_grid):
        try:
            SWEOS_eq = SpanWagnerEOS_EquilibriumPhase(None, P)
            rho_sat_gas[i] = SWEOS_eq.rho_sat_gas
        except Exception:
           rho_sat_gas[i] = np.nan

    np.savez(out_file, P_grid=P_grid, rho_sat_gas=rho_sat_gas)
    print("saved to:", out_file)
    return P_grid, rho_sat_gas

build_n2o_rho_sat_liq_table()
build_n2o_rho_sat_gas_table()