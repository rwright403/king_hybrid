from dataclasses import dataclass
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel
from src.utils.enum import FillType



@dataclass
class SimInputs:
    globals_kwargs: dict
    C: object
    ox_tank_kwargs: dict
    ox_inj_kwargs: dict
    fuel_tank_kwargs: dict | None
    cc_kwargs: dict
    nozzle_kwargs: dict


def build_kwargs(cfg):
    """Takes a config module and returns kwargs dicts depending on hybrid vs liquid."""

    # ------------------------
    # Global + environment
    # ------------------------
    globals_kwargs = dict(
        thrust_curve_graphs=cfg.thrust_curve_graphs,
        mode=cfg.mode,
        save_path=cfg.save_path,
        timestep=cfg.timestep,
        sim_time=cfg.sim_time,
        P_atm=cfg.P_atm,
        T_atm=cfg.T_atm,
        rho_atm=cfg.rho_atm,
    )

    models_kwargs = dict(
        ox_tank_model=cfg.ox_tank_model,
        fuel_tank_model=getattr(cfg, "fuel_tank_model", None),
        ox_inj_model=cfg.ox_inj_model,
        cc_model=cfg.cc_model,
        nozzle_model=cfg.nozzle_model,
    )

    # ------------------------
    # Build CEA object
    # ------------------------
    if cfg.fuel_properties is not None:  # hybrid case (custom solid fuel)
        add_new_fuel(cfg.fuel_name, cfg.fuel_properties)

    C = CEA_Obj(
        oxName=cfg.oxidizer_name,
        fuelName=cfg.fuel_name,
        pressure_units="Pa",
        isp_units="sec",
        cstar_units="m/s",
        temperature_units="K",
        sonic_velocity_units="m/s",
        enthalpy_units="kJ/kg",
        density_units="kg/m^3",
        specific_heat_units="kJ/kg-K",
    )

    # ------------------------
    # Oxidizer tank
    # ------------------------
    ox_tank_kwargs = dict(
        m_ox=cfg.m_ox,
        P_tank=cfg.P_ox_tank,
        P_cc=cfg.P_cc,
        P_atm=cfg.P_atm,
        T_atm=cfg.T_atm,
        rho_atm=cfg.rho_atm,
        V_tank=cfg.V_tank,
        diam_out=cfg.diam_out,
        diam_in=cfg.diam_in,
        rho_wall=cfg.rho_wall,
        k_w=cfg.k_w,
        volume_err_tol=cfg.volume_err_tol,
        P_dot_err_tol=cfg.P_dot_err_tol,
    )

    # ------------------------
    # Oxidizer injector
    # ------------------------
    ox_inj_kwargs = dict(
        Cd=cfg.Cd_inj,
        A_inj=cfg.A_inj_ox,
    )

    # ------------------------
    # Decide: Hybrid vs Liquid
    # ------------------------
    fuel_tank_kwargs = None
    cc_kwargs = None

    if hasattr(cfg, "fuel_tank_model"):  # Liquid engine path
        fuel_tank_kwargs = dict(
            m_fuel=cfg.m_fuel,
            m_pres=cfg.m_pres,
            P_pres_tank=cfg.P_pres_tank,
            P_atm=cfg.P_atm,
            T_atm=cfg.T_atm,
            rho_atm=cfg.rho_atm,
            V_fuel_tank=cfg.V_tank,
            V_pres_tank=cfg.V_pres_tank,
            diam_out=getattr(cfg, "diam_out_fuel", cfg.diam_out),
            diam_in=getattr(cfg, "diam_in_fuel", cfg.diam_in),
            rho_wall=getattr(cfg, "rho_wall_fuel", cfg.rho_wall),
            k_w=getattr(cfg, "k_w_fuel", cfg.k_w),
            Cd=cfg.Cd_inj,
            A_inj=cfg.A_inj_fuel,
            fuel_str=cfg.fuel_str,
            pres_str=cfg.pres_str,
            fill_type=FillType.ISOTHERMAL,
        )

        cc_kwargs = dict(
            L_star=cfg.L_star,
            C=C,   # no regression terms for liquid chamber
        )

    else:  # Hybrid engine path
        cc_kwargs = dict(
            L_star=cfg.L_star,
            m_fuel_i=cfg.m_fuel_i,
            rho_fuel=cfg.rho_fuel,
            a=cfg.a_reg,
            n=cfg.n_reg,
            L=cfg.L_port,
            A_port_i=cfg.A_port_i,
            C=C,
        )

    # ------------------------
    # Nozzle
    # ------------------------
    nozzle_kwargs = dict(
        d_throat=cfg.d_throat,
        expratio=cfg.expratio,
        P_atm=cfg.P_atm,
        C=C,
    )

    return dict(
        globals_kwargs=globals_kwargs,
        models_kwargs=models_kwargs,
        C=C,
        ox_tank_kwargs=ox_tank_kwargs,
        ox_inj_kwargs=ox_inj_kwargs,
        fuel_tank_kwargs=fuel_tank_kwargs,
        cc_kwargs=cc_kwargs,
        nozzle_kwargs=nozzle_kwargs,
    )
