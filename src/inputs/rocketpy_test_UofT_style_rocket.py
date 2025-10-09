# ------------------------
# Std Engine Input File (Converted Test Case)
# ------------------------

"""
NOTE: THIS IS NOT A REAL ROCKET
TOOK U OF T ENGINE DATA FROM 2022
AND PUT IT IN THEIR 2024 ROCKET PUBLISHED ON ROCKETPY
THIS IS JUST TO TEST THE FLIGHT SIM/ENGINE SIM COUPLING
"""

# Models
ox_tank_model = 2   
ox_inj_model  = 3   
cc_model      = 1   # hybrid_cc_w_fuel_grain
nozzle_model  = 1  

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 0.005       # [s]
sim_time = 6.0       # [s]
P_atm    = 101325     # [Pa]
T_atm    = 284.0      # [K]
rho_atm  = 1.225      # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------
fuel_name       = "paraffin" # NOTE: had issues with the rocketcea add fuel method so I modified library and added this card
oxidizer_name   = "N2O"

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 8.2           # [kg]
P_ox_tank   = 3920207.656   # [Pa]
V_ox_tank      = 0.012         # [m^3]
diam_out    = 0.1397        # [m]
diam_in     = 0.128524      # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]
volume_err_tol = 0.001       # from all_error
P_dot_err_tol = None #TODO: RM

# ------------------------
# Injector parameters
# ------------------------
Cd_inj_ox   = 0.611 # UTAT theoretically determined cd published in their report
A_inj_ox = 6.54E-05   # [m^2]


# ------------------------
# Chamber (hybrid regression model)
# ------------------------
V_pre_post_cc       = 0.00118174778  # [m^2]
m_fuel_i            = 1.5          # [kg]
rho_fuel            = 900.0        # [kg/m^3]
a_reg               = 0.155e-3     # [m/s*(kg/s)^n]
n_reg               = 0.5
L_port              = 0.3852333    # [m]
A_port              = 0.0038319753 # [m^2]

# ------------------------
# Nozzle
# ------------------------
d_throat = 0.037084      # [m]
expratio = 4.894398574   # [m]

# ------------------------
# Validation Files
# ------------------------
validation_files={
        "P_cc": "src/inputs/UofT_Deliverance_II/UofT_Deliverance_II_CC_Pressure.csv",
        "P_ox_tank": "src/inputs/UofT_Deliverance_II/UofT_Deliverance_II_Tank_Pressure.csv",
        "thrust": "src/inputs/UofT_Deliverance_II/UofT_Deliverance_II_Thrust.csv"
}

A_inj_fuel = None          # hybrid has no fuel injector




# ==========================
# ROCKETPY FLIGHT DESVAR
# ==========================
import datetime
from rocketpy import Environment

# --- LAUNCH CANADA TIMMINS ONT. LAUNCH PAD ---
"""
From Rocketpy Delivarance Mk. IV example flight (Launched at LC 2024)
(no wind, velo components set to 0)
https://docs.rocketpy.org/en/latest/examples/defiance_flight_sim.html
"""
flight_date = datetime.date(2024, 8, 24)
env = Environment(latitude=47.966527, longitude=-81.87413, elevation=1383.4)
env.set_date((flight_date.year, flight_date.month, flight_date.day, 0))
env.set_atmospheric_model(type="custom_atmosphere", wind_v=0.0, wind_u=0.0)

rail_length = 10
inclination = 85
heading = 90

# --- Rocket Geometry and Mass Properties ---
fuselage_radius = 0.07
rkt_dry_mass = 37.211
rkt_dry_inertia = (94.14, 94.14, 0.09)
rkt_dry_cg = 3.29
power_off_drag = "src/inputs/UofT_Deliverance_II/UofT_Style_Defiance_DragCurve.csv"
power_on_drag = "src/inputs/UofT_Deliverance_II/UofT_Style_Defiance_DragCurve.csv"
rkt_csys = "tail_to_nose"

upper_launch_lug_pos = 1.25
lower_launch_lug_pos = 0.5
launch_lug_angular_pos = 45

nose_length = 0.563
nose_kind = "vonKarman"
nose_position = 4.947

fins_n = 3
fins_span = 0.115
fins_root_chord = 0.4
fins_tip_chord = 0.2
fins_position = 0.175

### NOTE: is this cg or location actually??
ox_tank_cg = 2.2
engine_cg = 0

# --- Chamber / Motor Properties ---
cc_cg = 0.78
cc_dry_inertia = (1.801, 1.801, 0.0305)
cc_dry_mass = 13.832
nozzle_pos = 0.0
