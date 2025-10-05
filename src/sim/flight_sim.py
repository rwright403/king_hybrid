import os
import pandas as pd

from rocketpy import Rocket, Flight, HybridMotor, LiquidMotor, MassFlowRateBasedTank


def flight_sim(kwargs):
    """
    Runs a RocketPy flight simulation using propellant data
    from the propulsion simulation results.
    """


    ### Setup rocket and aero surfaces in this scope, then address motor!
    rocket = Rocket(
        radius=kwargs["rocketpy_rocket_kwargs"]["fuse_radius"],
        mass=kwargs["rocketpy_rocket_kwargs"]["rkt_motorless_mass"],
        inertia=kwargs["rocketpy_rocket_kwargs"]["rkt_motorless_inertia"],
        center_of_mass_without_motor=kwargs["rocketpy_rocket_kwargs"]["rkt_motorless_cg"],
        power_off_drag=kwargs["rocketpy_rocket_kwargs"]["power_off_drag"],
        power_on_drag=kwargs["rocketpy_rocket_kwargs"]["power_on_drag"],
        coordinate_system_orientation=kwargs["rocketpy_rocket_kwargs"]["rkt_csys"],
    )

    rocket.add_nose(
        length=kwargs["rocketpy_rocket_kwargs"]["nose_length"],
        kind=kwargs["rocketpy_rocket_kwargs"]["nose_kind"], 
        position=kwargs["rocketpy_rocket_kwargs"]["nose_position"],
        )

    rocket.add_trapezoidal_fins(
        n=kwargs["rocketpy_rocket_kwargs"]["fins_n"], 
        span=kwargs["rocketpy_rocket_kwargs"]["fins_span"], 
        root_chord=kwargs["rocketpy_rocket_kwargs"]["fins_root_chord"], 
        tip_chord=kwargs["rocketpy_rocket_kwargs"]["fins_tip_chord"], 
        position=kwargs["rocketpy_rocket_kwargs"]["fins_position"]
        )

    engine = None
#NOTE: IS IT EMPTY OR NONE?
    if kwargs.get("fuel_tank_kwargs") is not None: #We are flying liquid

        fuel_tank = MassFlowRateBasedTank(**kwargs["rocketpy_fuel_tank_kwargs"])
        engine = LiquidMotor(**kwargs["rocketpy_cc_kwargs"])

        engine.add_tank(fuel_tank, position=kwargs["rocketpy_rocket_kwargs"]["fuel_tank_cg"])

    else: #we are flying hybrid
        engine = HybridMotor(**kwargs["rocketpy_cc_kwargs"])


    ox_tank = MassFlowRateBasedTank(**kwargs["rocketpy_ox_tank_kwargs"])
    engine.add_tank(ox_tank, position=kwargs["rocketpy_rocket_kwargs"]["ox_tank_cg"])

    rocket.add_motor(engine, position=kwargs["rocketpy_rocket_kwargs"]["engine_cg"])

    rocket.set_rail_buttons(
        upper_button_position=kwargs["rocketpy_rocket_kwargs"]["upper_launch_lug_pos"],
        lower_button_position=kwargs["rocketpy_rocket_kwargs"]["lower_launch_lug_pos"],
        angular_position=kwargs["rocketpy_rocket_kwargs"]["launch_lug_angular_pos"],
    )

    rocket.draw()

    flight = Flight(
        rocket=rocket,
        environment=kwargs["rocketpy_launchpad_kwargs"]["env"],
        rail_length=kwargs["rocketpy_launchpad_kwargs"]["rail_length"],
        inclination=kwargs["rocketpy_launchpad_kwargs"]["inclination"],
        heading=kwargs["rocketpy_launchpad_kwargs"]["heading"],
    )
    flight.prints.out_of_rail_conditions()
    flight.plots.trajectory_3d()

    return flight