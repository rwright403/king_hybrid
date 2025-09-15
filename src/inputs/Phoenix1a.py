# ------------------------
# Std Engine Input File (Converted Test Case - Hybrid Rocket)
# ------------------------

# Models
ox_tank_model = 1   # bens_ox_tank
ox_inj_model  = 1   # SPI injector
cc_model      = 1   # hybrid_cc_w_fuel_grain
nozzle_model  = 1   # nozzle model

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 0.05       # [s]
sim_time = 8.0        # [s]
P_atm    = 101325     # [Pa]
T_atm    = None       # not given
rho_atm  = None       # not given

# ------------------------
# Propellant properties
# ------------------------
fuel_name       = "paraffin"
fuel_properties = """fuel paraffin  C 20   H 42    wt%=100.00
h,KJ/mol=-13313.49  t(k)=298.15   rho,kg/m3=900
"""
oxidizer_name   = "N2O"

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 18.30872404   # [kg]
P_ox_tank   = 6.5e6         # [Pa]
P_cc        = None          # not given
V_tank      = 0.043         # [m^3]
diam_out    = None
diam_in     = None
rho_wall    = None
k_w         = None
volume_err_tol = 0.01
P_dot_err_tol  = None

# ------------------------
# Injector parameters
# ------------------------
Cd_inj   = None               # not given separately
A_inj_ox = 2.243097155e-5     # [m^2] effective Cd*A provided
A_inj_fuel = None             # hybrid has no fuel injector

# ------------------------
# Chamber (hybrid regression model)
# ------------------------
L_star    = None
m_fuel_i  = 6.173983      # [kg]
rho_fuel  = 900.0         # [kg/m^3]
a_reg     = 0.155e-3      # [m/s*(kg/s)^n]
n_reg     = 0.5
L_port    = 0.4           # [m]
A_port  = 0.0019634954  # [m^2]

# ------------------------
# Nozzle
# ------------------------
# Provided areas, leave as None unless we compute
d_throat = None
expratio = None

# ------------------------
# Experimental validation data
# ------------------------
exp_thrust_file_path   = r'./src/inputs/bens_validation_data/Phoenix_1a/Phoenix_1a_Thrust.csv'
exp_p_cc_file_path     = r'./src/inputs/bens_validation_data/Phoenix_1a/Phoenix_1a_CC_Pressure.csv'
exp_p_ox_tank_file_path = r'./src/inputs/bens_validation_data/Phoenix_1a/Phoenix_1a_Tank_Pressure.csv'
