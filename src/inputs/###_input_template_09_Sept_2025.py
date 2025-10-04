# ************************
# INPUT FILE TEMPLATE
# ************************


# ========================
# ENGINE MODEL DESVAR
# ========================

# Models
ox_tank_model = 2   #  2 = non-equilibrium, 1 = equilibrium (numerically unstable, not recommended)
fuel_tank_model = 1 #NOTE: ONLY LRE
ox_inj_model  = 1   # injector model selection
cc_model      = 1   # chamber model
nozzle_model  = 1   # nozzle model


# Global settings
thrust_curve_graphs = True
mode = "full_stack"

# ------------------------
# Simulation control
# ------------------------
# Options:
#   "ox_tank"   → simulate oxidizer tank only
#   "fuel_tank" → simulate fuel/pressurant tank only
#   "full_stack" → simulate full engine (default)
# NOTE: must run full_stack for flight sim analysis
mode = "full_stack"
save_path = None #TODO: RM

# Engine Environment
TIMESTEP = 1e-2       # [s]
sim_time = 3.0        # [s]
P_atm    = 101325     # [Pa]
T_atm    = 298.15     # [K]
rho_atm  = 1.225      # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------ #TODO: MAKE BETTER
fuel_name       = "HTPB"
fuel_properties = "HTPB C=0.9 H=0.1"
oxidizer_name   = "N2O"
fuel_str        = "Ethanol"   # for liquid-fuel config
pres_str        = "Nitrogen"  # pressurant gas (alt: "Helium")

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 1         # [kg]
P_ox_tank   = 1      # [Pa]
V_ox_tank   = 1       # [m^3]
diam_out    = 1      # [m]
diam_in     = 1        # [m]
rho_wall    = 1          # [kg/m^3]
k_w         = 1           # [W/(m K)]
volume_err_tol = 1       # [m^3]
P_dot_err_tol  = None

# ------------------------
# Ox Injector parameters
# ------------------------
Cd_inj_ox   = 1
A_inj_ox = 1     # [m^2]


# ------------------------
# Fuel tank & pressurant
# ------------------------ #NOTE: LRE ONLY
m_fuel       = 1            # [kg]
m_pres       = 1            # [kg] pressurant mass (estimated)
P_pres_tank  = 1           # [Pa] (from P_fueltank)
V_fuel_tank  = 1          # [m^3] fuel tank volume
V_pres_tank  = 1        # [m^3]  #NOTE: NO EXT PRESSURANT TANK!
diam_out_fuel = 1          # [m]
diam_in_fuel  = 1           # [m]
rho_wall    = 1          # [kg/m^3]
k_w         = 1            # [W/(m K)]


# ------------------------
# Fuel Injector parameters
# ------------------------ #NOTE: LRE ONLY
Cd_inj_fuel = 1         # [-]
A_inj_fuel  = 1         # [m^2]

# ------------------------
# Liquid Rocket Engine Chamber
# ------------------------
V_cc = 1             # [m^3]


# ------------------------
# Chamber (fuel regression)
# ------------------------
L_star     = 1         # [m]
m_fuel_i   = 1        # [kg] initial solid fuel
rho_fuel   = 1      # [kg/m^3]
a_reg      = 1     # [m/s*(kg/s)^n] regression rate coeff
n_reg      = 1        # [-] regression rate exponent
L_port     = 1        # [m] port length
A_port      = 1      # [m^2] initial port cross-sectional area

# ------------------------
# Nozzle
# ------------------------
d_throat   = 1     # [m]
expratio   = 1        # A_exit / A_throat


# ------------------------
# Validation Files (for engine model)
# ------------------------
validation_files={
    "P_cc": "src/inputs/...",
    "P_ox_tank": "src/inputs/...",
    "P_fuel_tank": "src/inputs/...", # NOTE: LRE ONLY
    "thrust": "src/inputs/...",
}



# ==========================
# ROCKETPY FLIGHT DESVAR
# ==========================
from datetime import datetime
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

rail_length = 1
inclination = 1
heading = 1

# --- Rocket Geometry and Mass Properties ---
fuselage_radius = 1
rkt_dry_mass = 1
rkt_dry_inertia = 1
rkt_dry_cg = 1
power_off_drag = 1
power_on_drag = 1
rkt_csys = 1

upper_launch_lug_pos = 1
lower_launch_lug_pos = 1
launch_lug_angular_pos = 1

nose_length = 1
nose_kind = 1
nose_position = 1

fins_n = 1
fins_span = 1
fins_root_chord = 1
fins_tip_chord = 1
fins_position = 1

### NOTE: is this cg or location actually??
ox_tank_cg = 1
fuel_tank_cg = 1
engine_cg = 1

# --- Chamber / Motor Properties ---
cc_cg = 1
cc_dry_inertia = 1
cc_dry_mass = 1
noz_pos = 1
