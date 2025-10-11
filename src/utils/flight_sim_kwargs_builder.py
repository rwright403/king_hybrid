import os
import numpy as np
import CoolProp.CoolProp as CP
from rocketpy import Fluid, CylindricalTank
from src.utils.model_registry import *

def get_model(kind: str, code: int):
    if kind == "D":
        module_path = DRAG_MODEL_MAP.get(code)
    else:
        raise ValueError(f"Unknown model kind {kind}")
    
    if module_path is None:
        raise ValueError(f"Model {code} is not defined in {kind} map")

    module_path, class_name = module_path.rsplit(".", 1)
    module = importlib.import_module(module_path)
    return getattr(module, class_name)


def build_flight_sim_kwargs(input_file, cfg):
    """Takes a config module and returns kwargs dicts depending on hybrid vs liquid."""

    # --- Construct dynamic path to engine output folder ---
    engine_output_dir = os.path.join("src", "results", str(input_file))

    ENGINE_MODEL_OUTPUT_MAP = {
        "m_dot_ox": os.path.join(engine_output_dir, "m_dot_ox.csv"),
        "m_dot_fuel": os.path.join(engine_output_dir, "m_dot_fuel.csv"),
        "thrust": os.path.join(engine_output_dir, "thrust.csv"),
    }

    rocketpy_fuel_tank_kwargs = None

    # ------------------------
    # Environment
    # ------------------------
    rocketpy_launchpad_kwargs = dict(
        env=cfg.env,
        rail_length=cfg.rail_length,
        inclination=cfg.inclination,
        heading=cfg.heading
    ) 


    # ------------------------
    # Drag Model Setup
    # ------------------------
    #NOTE: need to build drag model in this scope

    Cd_power_off = cfg.power_off_drag
    Cd_power_on = cfg.power_on_drag       

    if cfg.drag_model != 1:

        dragClass = get_model("D", cfg.drag_model)
            
        if cfg.drag_model == 2:

            """
            drag_kwargs = dict(
                L = cfg.L,                      
                D = cfg.D,                      
                fin_area = cfg.fin_area,
                ref_area = cfg.ref_area,
                N = cfg.N_fins,                 
                rL = cfg.r_le,                  
                GammaL = cfg.sweep_angle,       
                Cr = cfg.c_root,  
                hr = cfg.fin_te_thickness, 
                tr = cfg.fin_thickness,  
                S = cfg.fin_semispan,   
                Ar_fins = cfg.fin_ref_area, 
                R_s = cfg.surface_roughness,
                rho = cfg.rho_atm,
            )
            """
            drag_kwargs = dict(
                x=2
            )

        drag_model = dragClass()

        Cd_power_off = drag_model.make_drag_model(drag_kwargs)
        Cd_power_on = drag_model.make_drag_model(drag_kwargs)

    


    # ------------------------
    # Rocket
    # ------------------------ 
    rocketpy_rocket_kwargs = dict(
        fuse_radius=cfg.fuselage_radius,
        rkt_motorless_mass=cfg.rkt_dry_mass,
        rkt_motorless_inertia=cfg.rkt_dry_inertia,
        rkt_motorless_cg=cfg.rkt_dry_cg,
        power_off_drag=Cd_power_off,
        power_on_drag=Cd_power_on,
        rkt_csys="tail_to_nose",

        upper_launch_lug_pos=cfg.upper_launch_lug_pos,
        lower_launch_lug_pos=cfg.lower_launch_lug_pos,
        launch_lug_angular_pos=cfg.launch_lug_angular_pos,

        nose_length=cfg.nose_length, 
        nose_kind=cfg.nose_kind, 
        nose_position=cfg.nose_position,

        fins_n=cfg.fins_n, 
        fins_span=cfg.fins_span, 
        fins_root_chord=cfg.fins_root_chord, 
        fins_tip_chord=cfg.fins_tip_chord, 
        fins_position=cfg.fins_position,

        ox_tank_cg=cfg.ox_tank_cg,

        engine_cg=cfg.engine_cg

    )
    #NOTE: WE ADD "fuel_tank_cg" if there is a fuel tank below with the rocketpy add tank function



    # ------------------------
    # Oxidizer tank
    # ------------------------ 
    #    
    """
    NOTE: from rocketpy documentation for MassFlowRateBasedTank:
    If a .csv file is given, it must have two columns, the first one being time in seconds and the second one being the mass flow rate in kg/s.
    """
    m_dot_ox_filepath = ENGINE_MODEL_OUTPUT_MAP["m_dot_ox"]
        
    ox_tank_height = (cfg.V_ox_tank)/(0.25*np.pi*(cfg.diam_in**2))

    rho_ox_liq = CP.PropsSI('D', 'P', cfg.P_ox_tank, 'Q', 0, "N2O")
    rho_ox_gas = CP.PropsSI('D', 'P', cfg.P_ox_tank, 'Q', 1, "N2O")
    rho_bulk_tank = cfg.m_ox/cfg.V_ox_tank


    x_tank = ( (1/rho_bulk_tank)-(1/rho_ox_liq) ) / ( (1/rho_ox_gas)-(1/rho_ox_liq) )
    m_ox_gas = x_tank*cfg.m_ox
    m_ox_liq = cfg.m_ox-m_ox_gas


    # RocketPy represents quantities like volume and height as Function objects (not plain floats).
    # When initializing the tank, it composes the geometry inverse function h(V) with the current
    # fluid volume Function V(t). Before doing so, RocketPy verifies that the *entire range* of V(t)
    # (its min and max possible values) lies within the geometry’s valid volume domain.
    #
    # In our case, the total fluid volume computed from m_liq/ρ_liq + m_gas/ρ_gas is just slightly
    # larger (≈1e-8 m³) than the tank’s geometric volume due to floating-point rounding.
    # That makes RocketPy think V(t) ∈ [0.011994, 0.012], while the geometry only supports
    # [0.0, 0.0119999], causing a domain violation.
    #
    # Simple fix: slightly perturb the tank height (e.g., +1%) so the geometric volume is
    # marginally larger than the computed fluid volume. This keeps the fluid Function’s image
    # strictly inside the geometry’s domain and prevents RocketPy’s inverse-volume error.
    ox_tank_height*=1.01


    rocketpy_ox_tank_kwargs = dict(
            name="ox_tank",
            geometry=CylindricalTank( (0.5*cfg.diam_in), ox_tank_height), #radius, height
            flux_time=cfg.sim_time,
            liquid=Fluid(name="ox_liq", density=rho_ox_liq),
            gas=Fluid(name="ox_gas", density=rho_ox_gas),
            initial_liquid_mass=m_ox_liq,
            initial_gas_mass=m_ox_gas,
            liquid_mass_flow_rate_in=0.0,
            liquid_mass_flow_rate_out=m_dot_ox_filepath,
            gas_mass_flow_rate_in=0.0,
            gas_mass_flow_rate_out=0.0,
        )



    # ------------------------
    # Decide: Hybrid vs Liquid
    # ------------------------

    thrust_filepath = ENGINE_MODEL_OUTPUT_MAP["thrust"]

    A_exit = cfg.expratio*0.25*np.pi*(cfg.d_throat**2)
    r_nozzle_exit = np.sqrt(A_exit/np.pi)

    fuel_tank_model=getattr(cfg, "fuel_tank_model", None)
    if fuel_tank_model is not None:  # Liquid engine path

        m_dot_fuel_filepath = ENGINE_MODEL_OUTPUT_MAP["m_dot_fuel"]
        
        fuel_tank_height = (cfg.V_fuel_tank+cfg.V_pres_tank)/(0.25*np.pi*(cfg.diam_in**2))

        rho_fuel = CP.PropsSI('D', 'P', cfg.P_pres_tank, 'T', cfg.T_atm, cfg.fuel_str)
        rho_pres = CP.PropsSI('D', 'P', cfg.P_pres_tank, 'T', cfg.T_atm, cfg.pres_str)

        rocketpy_fuel_tank_kwargs = dict(
            name="fuel_tank",
            geometry=CylindricalTank( (0.5*cfg.diam_in_fuel), fuel_tank_height), #radius, height
            flux_time=cfg.sim_time,
            liquid=Fluid(name="fuel_liq", density=rho_fuel),
            gas=Fluid(name="pres_gas", density=rho_pres),
            initial_liquid_mass=cfg.m_fuel,
            initial_gas_mass=cfg.m_gas,
            liquid_mass_flow_rate_in=0.0,
            liquid_mass_flow_rate_out=m_dot_fuel_filepath,
            gas_mass_flow_rate_in=0.0,
            gas_mass_flow_rate_out=0.0,
        )

        rocketpy_rocket_kwargs["fuel_tank_cg"] = cfg.fuel_tank_cg

        rocketpy_cc_kwargs = dict(
            thrust_source=thrust_filepath,
            center_of_dry_mass_position=cfg.cc_cg,
            dry_inertia=cfg.cc_dry_inertia,
            dry_mass=cfg.cc_dry_mass,
            burn_time=cfg.sim_time,
            nozzle_radius=r_nozzle_exit, #NOTE: Radius of motor nozzle outlet in meters.
            nozzle_position=cfg.nozzle_pos,
            coordinate_system_orientation="nozzle_to_combustion_chamber",
        )



    else: #its a hybrid!

        A_fuel_grain_od = (cfg.m_fuel_i / (cfg.rho_fuel * cfg.L_port)) + cfg.A_port
        r_fuel_grain_outer = np.sqrt(A_fuel_grain_od/np.pi)
        r_fuel_grain_inner = np.sqrt(cfg.A_port/np.pi)


        rocketpy_cc_kwargs = dict(
            thrust_source=thrust_filepath,
            dry_mass=cfg.cc_dry_mass,
            dry_inertia=cfg.cc_dry_inertia,
            center_of_dry_mass_position=cfg.cc_cg,
            reshape_thrust_curve=False,
            grain_number=1,
            grain_separation=0,
            grain_outer_radius=r_fuel_grain_outer,
            grain_initial_inner_radius=r_fuel_grain_inner,
            grain_initial_height=cfg.L_port,
            grain_density=cfg.rho_fuel,
            nozzle_radius=r_nozzle_exit,
            throat_radius=(0.5*cfg.d_throat),
            interpolation_method="linear",
            grains_center_of_mass_position=cfg.cc_cg, #NOTE: ASSUME CG GRAIN = CG CHAMBER
            coordinate_system_orientation="nozzle_to_combustion_chamber",
        )



    return {
        "rocketpy_launchpad_kwargs": rocketpy_launchpad_kwargs,
        "rocketpy_rocket_kwargs": rocketpy_rocket_kwargs,
        "rocketpy_ox_tank_kwargs": rocketpy_ox_tank_kwargs,
        "rocketpy_fuel_tank_kwargs": rocketpy_fuel_tank_kwargs,
        "rocketpy_cc_kwargs": rocketpy_cc_kwargs,
    }
