from dataclasses import dataclass
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel
from src.utils.enum import FillType



def build_kwargs(cfg):
    """Takes a config module and returns kwargs dicts depending on hybrid vs liquid."""

    # ------------------------
    # Global + environment
    # ------------------------
    globals_kwargs = dict(
        thrust_curve_graphs=getattr(cfg, "thrust_curve_graphs", None),
        mode=getattr(cfg, "mode", None),
        save_path=getattr(cfg, "save_path", None),
        timestep=getattr(cfg, "timestep", None),
        sim_time=getattr(cfg, "sim_time", None),
        P_atm=getattr(cfg, "P_atm", None),
        T_atm=getattr(cfg, "T_atm", None),
        rho_atm=getattr(cfg, "rho_atm", None),
    )


    models_kwargs = dict(
        ox_tank_model=getattr(cfg, "ox_tank_model", None),
        fuel_tank_model=getattr(cfg, "fuel_tank_model", None),
        ox_inj_model=getattr(cfg, "ox_inj_model", None),
        cc_model=getattr(cfg, "cc_model", None),
        nozzle_model=getattr(cfg, "nozzle_model", None),
    )

    # ------------------------
    # Build CEA object (only if CC/nozzle exists)
    # ------------------------
    C = None
    if (models_kwargs["cc_model"] or models_kwargs["nozzle_model"]) \
   and getattr(cfg, "fuel_name", None) is not None \
   and getattr(cfg, "oxidizer_name", None) is not None:

        """
        if getattr(cfg, "fuel_properties", None) is not None:
            add_new_fuel(cfg.fuel_name, cfg.fuel_properties)
        """
            
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

        #print("kwargs builder MW: unit check: ", C.get_Chamber_MolWt_gamma(1e6, 6.0, 40)[0]) #: return the tuple (mw, gam)

    # ------------------------
    # Oxidizer tank
    # ------------------------
    ox_tank_kwargs = None
    if models_kwargs["ox_tank_model"] == 1:
        ox_tank_kwargs = dict(
            m_ox=cfg.m_ox,
            V_tank=cfg.V_tank,
            P_tank=cfg.P_ox_tank,
            P_atm=cfg.P_atm,
            all_error=cfg.volume_err_tol,
        )
    elif models_kwargs["ox_tank_model"] == 2:
        ox_tank_kwargs = dict(
            m_ox=cfg.m_ox,
            P_tank=cfg.P_ox_tank,
            P_atm=cfg.P_atm,
            T_atm=cfg.T_atm,
            rho_atm=cfg.rho_atm,
            V_tank=cfg.V_ox_tank,
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
    ox_inj_kwargs = None
    if models_kwargs["ox_inj_model"] is not None:
        ox_inj_kwargs = dict(
            Cd=cfg.Cd_inj_ox,
            A_inj=cfg.A_inj_ox,
        )

    # ------------------------
    # Decide: Hybrid vs Liquid
    # ------------------------
    fuel_tank_kwargs = None
    cc_kwargs = None


    if models_kwargs["fuel_tank_model"] is not None:  # Liquid engine path
        fuel_tank_kwargs = dict(
            m_fuel=cfg.m_fuel,
            m_pres=cfg.m_pres,
            P_pres_tank=cfg.P_pres_tank,
            P_atm=cfg.P_atm,
            T_atm=cfg.T_atm,
            rho_atm=cfg.rho_atm,
            V_fuel_tank=cfg.V_fuel_tank,
            V_pres_tank=cfg.V_pres_tank,
            diam_out=getattr(cfg, "diam_out_fuel", cfg.diam_out),
            diam_in=getattr(cfg, "diam_in_fuel", cfg.diam_in),
            rho_wall=getattr(cfg, "rho_wall_fuel", cfg.rho_wall),
            k_w=getattr(cfg, "k_w_fuel", cfg.k_w),
            Cd=cfg.Cd_inj_fuel,
            A_inj=cfg.A_inj_fuel,
            fuel_str=cfg.fuel_str,
            pres_str=cfg.pres_str,
            fill_type=FillType.ISOTHERMAL,
        )

    if models_kwargs["cc_model"] == 1:  # Hybrid engine path
        cc_kwargs = dict(
            V_pre_post_cc=cfg.V_pre_post_cc,
            m_fuel_i=cfg.m_fuel_i,
            rho_fuel=cfg.rho_fuel,
            a=cfg.a_reg,
            n=cfg.n_reg,
            L=cfg.L_port,
            A_port=cfg.A_port,
            C=C,
        )
    elif models_kwargs["cc_model"] == 2:  # Liquid Engine Path
        cc_kwargs = dict(
            V_cc=cfg.V_cc,
            C=C,
        )
    elif models_kwargs["cc_model"] == 3:
        cc_kwargs = dict(
            filepath=cfg.validation_files["P_cc"],
            timestep=cfg.timestep
        )

    # ------------------------
    # Nozzle
    # ------------------------
    nozzle_kwargs = None
    if models_kwargs["nozzle_model"] == 1:
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
