from dataclasses import dataclass
import numpy as np
import CoolProp.CoolProp as CP
from src.utils.build_rocketpy_input_csv import build_rocketpy_input_csv
from rocketpy import Fluid, CylindricalTank
from src.utils.build_rocketpy_input_csv import build_rocketpy_input_csv

"""
@dataclass
class FlightSimRocketInputs:
    rocketpy_ox_tank_kwargs: dict | None
    rocketpy_fuel_tank_kwargs: dict | None
    rocketpy_cc_kwargs: dict | None
"""

def build_flight_sim_kwargs(cfg, prop_results):
    """Takes a config module and returns kwargs dicts depending on hybrid vs liquid."""

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
    # Rocket
    # ------------------------ 
    rocketpy_rocket_kwargs = dict(
        fuse_radius=cfg.fuselage_radius,
        rkt_motorless_mass=cfg.rkt_dry_mass,
        rkt_motorless_inertia=cfg.rkt_dry_inertia,
        rkt_motorless_cg=cfg.rkt_dry_cg,
        power_off_drag=cfg.power_off_drag,
        power_on_drag=cfg.power_on_drag,
        rkt_csys="tail_to_nose",

        upper_launch_lug_pos=cfg.upper_launch_lug_pos,
        lower_launch_lug_pos=cfg.lower_launch_lug_pos,
        launch_lug_angular_pos=cfg.launch_lug_angular_pos,

        nose_length=cfg.nose_length, 
        nose_kind=cfg.nose_kind, 
        nose_position=cfg.nose_pos,

        fins_n=cfg.fins_n, 
        fins_span=cfg.fins_span, 
        fins_root_chord=cfg.fins_root_chord, 
        fins_tip_chord=cfg.fins_tip_chord, 
        fins_position=cfg.fins_position,

        ox_tank_cg=cfg.ox_tank_cg,

        engine_cg=cfg.engine_cg

    )
    #NOTE: WE ADD "fuel_tank_cg" if there is a fuel tank



    # ------------------------
    # Oxidizer tank
    # ------------------------ 
    #    
    """
    NOTE: from rocketpy documentation for MassFlowRateBasedTank:
    If a .csv file is given, it must have two columns, the first one being time in seconds and the second one being the mass flow rate in kg/s.
    """
    m_dot_ox_filepath = build_rocketpy_input_csv(prop_results, "m_dot_ox")
        
    ox_tank_height = (cfg.V_ox_tank)/(0.25*np.pi*(cfg.diam_in**2))

    rho_ox_liq = CP.PropsSI('D', 'P', cfg.P_ox_tank, 'T', cfg.T_atm, "N2O")
    rho_ox_gas = CP.PropsSI('D', 'P', cfg.P_ox_tank, 'T', cfg.T_atm, "N2O")
    rho_bulk_tank = cfg.m_ox/cfg.V_ox_tank
    x_tank = ( (1/rho_bulk_tank)-(1/rho_ox_liq) ) / ( (1/rho_ox_gas)-(1/rho_ox_liq) )
    m_ox_gas = x_tank*cfg.m_ox
    m_ox_liq = cfg.m_ox-m_ox_gas


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

    thrust_filepath = build_rocketpy_input_csv(prop_results, "thrust")

    A_exit = cfg.expratio*0.25*np.pi*(cfg.d_throat**2)
    r_noz_exit = np.sqrt(A_exit/np.pi)

    fuel_tank_model=getattr(cfg, "fuel_tank_model", None)
    if fuel_tank_model is not None:  # Liquid engine path

        m_dot_fuel_filepath = build_rocketpy_input_csv(prop_results, "m_dot_fuel")
        
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
            nozzle_radius=r_noz_exit, #NOTE: Radius of motor nozzle outlet in meters.
            nozzle_position=cfg.noz_pos,
            coordinate_system_orientation="nozzle_to_combustion_chamber",
        )



    else: #its a hybrid!

        A_fuel_grain_od = (cfg.m_fuel_i / (cfg.rho_fuel * cfg.L)) + cfg.A_port
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
            grain_initial_height=cfg.L,
            grain_density=cfg.rho_fuel,
            nozzle_radius=r_noz_exit,
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
