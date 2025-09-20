"""
SETTING UP A SIM WITH ROCKETPY
DEFINE ENVIRONMENT --> DEFINE MOTOR --> DEFINE ROCKET
RUNNING SIM BY DEFINING FLIGHT OBJECT
PLOTTING AT THE END
"""

from src import constants
import numpy as np

from rocketpy import (
    Environment,
    Rocket,
    Flight,
    HybridMotor,
    Fluid,
    CylindricalTank,
    MassFlowRateBasedTank,

)

#1) def environment
env = Environment(latitude=constants.latitude, longitude=constants.latitude, elevation=constants.elevation) 

env.set_date(
    (constants.year, constants.month, constants.date, constants.hour)
)  # Hour given in UTC time #Arbitrary Date

env.set_atmospheric_model(type="standard_atmosphere", file="GFS") #THIS USES GFS ATMOSPHERIC MODEL, OTHER MODELS EXIST CAN BE LOOKED AT LATER
#env.info() #PRINTS ENVIRONMENT DATA

#2) Define ENGINE! START W TANK THEN CC

# Define the fluids
oxidizer_liq = Fluid(name="N2O_l", density=constants.rho_ox_liq)
oxidizer_gas = Fluid(name="N2O_g", density=constants.rho_ox_gas)

# Define tank geometry
tank_shape = CylindricalTank(constants.r_tank, constants.height_tank, spherical_caps=False)

# Define tank
oxidizer_tank = MassFlowRateBasedTank(
    name="oxidizer tank",
    geometry=tank_shape,
    flux_time=5.2,
    initial_liquid_mass=constants.m_ox*constants.fill_level,
    initial_gas_mass=constants.m_ox*(1-constants.fill_level),
    liquid_mass_flow_rate_in=0,
    liquid_mass_flow_rate_out=r'./src/m_dot_ox.csv',
    gas_mass_flow_rate_in=0,
    gas_mass_flow_rate_out=0,
    liquid=oxidizer_liq,
    gas=oxidizer_gas,
)

#COMBUSTION CHAMBER
example_hybrid = HybridMotor(
    thrust_source=r'./src/thrust.csv',
    dry_mass=2,
    dry_inertia=(0.125, 0.125, 0.002),
    nozzle_radius=63.36 / 2000,
    grain_number=4,
    grain_separation=0,
    grain_outer_radius=0.0723554901/2, #TODO: MAKE THIS A CONSTANT
    grain_initial_inner_radius=np.sqrt(constants.A_port/np.pi),
    grain_initial_height=constants.L,
    grain_density=constants.rho_fuel,
    grains_center_of_mass_position=constants.L/2,
    center_of_dry_mass_position=0.284,
    nozzle_position=0,
    burn_time=7,
    throat_radius=np.sqrt(constants.A_throat/np.pi),
)

#ADD THE OX TANK TO THE CC
example_hybrid.add_tank(
  tank = oxidizer_tank, position = 1.50115
)

#NOW WE HAVE OUR HYBRID, CAN SEE RESULTS W example_hybrid.all_info()

#NOW WE DEFINE THE ROCKET:
"""
Define the rocket itself by passing in the rockets dry mass, inertia, drag coefficient and radius;

Add a motor;

Add, if desired, aerodynamic surfaces;

Add, if desired, parachutes;

Set, if desired, rail guides;

See results.
"""
#CREATE ROCKET
XENIA2 = Rocket(
    radius=constants.rocket_fuselage_rad,
    mass=constants.rocket_dry_mass,
    inertia=(6.321, 6.321, 0.034),
    power_off_drag=r'./src/flight_sim/sample_drag.csv',
    power_on_drag=r'./src/flight_sim/sample_drag.csv',
    center_of_mass_without_motor=0,
    coordinate_system_orientation="tail_to_nose",
)


#NEED TO ADD THE HYBRID TO THE  XENIA 2 ROCKET
XENIA2.add_motor(example_hybrid, position=0)

"""
#can add rail guides to the rocket to get a better idea of off the rail velocity and stability!!!

rail_buttons = calisto.set_rail_buttons(
    upper_button_position=0.0818,
    lower_button_position=-0.6182,
    angular_position=45,
)
"""

#CAN ADD AERO SURFACES TO ROCKET
nose_cone = XENIA2.add_nose(
    length=0.55829, kind="von karman", position=3
)

fin_set = XENIA2.add_trapezoidal_fins(
    n=3,
    root_chord=0.120,
    tip_chord=0.060,
    span=0.110,
    position=-1.504956,
    cant_angle=0.5,
    airfoil=(r'./src/flight_sim/NACA0012-radians.csv',"radians"),
)

tail = XENIA2.add_tail(
    top_radius=0.0635, bottom_radius=0.0435, length=0.060, position=0
)

"""
#CAN ADD PARACHUTES!!!!!
main = XENIA2.add_parachute(
    name="main",
    cd_s=10.0,
    trigger=800,      # ejection altitude in meters
    sampling_rate=105,
    lag=1.5,
    noise=(0, 8.3, 0.5),
)

drogue = XENIA2.add_parachute(
    name="drogue",
    cd_s=1.0,
    trigger="apogee",  # ejection at apogee
    sampling_rate=105,
    lag=1.5,
    noise=(0, 8.3, 0.5),
)
"""

####SUGGESTED TO PLOT STATIC MARGIN TO SEE STABILITY. SIM WILL FAIL IF OVERSTABLE OR NEGATIVE
XENIA2.plots.static_margin()

####LASTLY TO RUN SIM WE NEED TO CREATE A FLIGHT OBJECT
test_flight = Flight(
    rocket=XENIA2, environment=env, rail_length=constants.launch_rail_length, inclination=constants.inclination, heading=constants.heading
    )

test_flight.all_info() #BIG DATA