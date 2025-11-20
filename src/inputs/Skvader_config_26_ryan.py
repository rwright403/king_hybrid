# ------------------------
# Std Engine Input File (Converted Test Case)
# ------------------------

import matplotlib
matplotlib.use("TkAgg")  # or "Qt5Agg" if you have Qt installed
import matplotlib.pyplot as plt

# Models
ox_tank_model = 2  
ox_inj_model  = 3  
cc_model      = 1   # hybrid_cc_w_fuel_grain
nozzle_model  = 1  
drag_model    = 1

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 0.005     # [s]
sim_time = 12      # [s]
P_atm    = 101325     # [Pa]
T_atm    = 273.15 + 22      # [K]
rho_atm  = 1.225      # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------
fuel_name       = "paraffin" # NOTE: had issues with the rocketcea add fuel method so I modified library and added this card
oxidizer_name   = "N2O"

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 17.7       # [kg]
P_ox_tank   = 53e5          # [Pa]


import CoolProp.CoolProp as CP
P_sat = CP.PropsSI('P', 'T', T_atm, 'Q', 0, 'N2O')
print(P_sat, P_ox_tank)

V_ox_tank   = 0.0315 #0.04 #0.0315      # [m^3] 0.75 Fill Fraction
diam_out    = 0.14          # [m]
diam_in     = 0.13365       # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]
volume_err_tol = 0.001      # from all_error
P_dot_err_tol = None #TODO: RM

# ------------------------
# Injector parameters
# ------------------------
# using [6]
Cd_inj_ox   = 0.57
n = 42
A_inj_ox = 0.25 * 3.14159 * n * (1.5e-3)**2   # [m^2]


# ------------------------
# Chamber (hybrid regression model) 
# ------------------------
V_pre_post_cc       = 0.00118174778  # [m^2]
m_fuel_i            = 2.36           # [kg]
rho_fuel            = 900.0          # [kg/m^3]
a_reg               = 0.155e-3       # [m/s*(kg/s)^n]
n_reg               = 0.5
L_port              = 0.5            # [m]
A_port              = 0.001864       # [m^2]

import numpy as np

V_total = m_fuel_i / rho_fuel
V_inner = A_port * L_port
A_outer = (V_total + V_inner) / L_port

D_outer = 2 * np.sqrt(A_outer / np.pi)
print("Outer diameter [in]:", 39.3701*D_outer)


# ------------------------
# Nozzle
# ------------------------
d_throat = 0.036673386    # [m]
expratio = 8.413859791     # [m]


# ------------------------
# Sizing wizard estimates (kept for reference, not used in runtime)
# ------------------------
min_TW_ratio = 11
Cd_est = 0.6
mass_fraction_estimate = 0.2657741984
characteristic_len = 1.0        #[m]
elevation = 225 # TODO: CHECK ALT MALAHAT RANGE VS TIMMINS
apogee_height = 3048*3          #[m]


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