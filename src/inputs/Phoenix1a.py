# ------------------------
# Std Engine Input File (Converted Test Case - Hybrid Rocket)
# ------------------------

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
timestep = 0.05       # [s]
sim_time = 8.0        # [s]
P_atm    = 101325     # [Pa]
T_atm    = 284.0      # [K]
rho_atm  = 1.225      # [kg/m^3]

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
Cd_inj   = 0.6               # not given separately
A_inj_ox = 2.243097155e-5     # [m^2] effective Cd*A provided
A_inj_fuel = None             # hybrid has no fuel injector

# ------------------------
# Chamber (hybrid regression model)
# ------------------------
V_pre_post_cc       = 0.00118174778  # [m^2]
m_fuel_i            = 6.173983      # [kg]
rho_fuel            = 900.0         # [kg/m^3]
a_reg               = 0.155e-3      # [m/s*(kg/s)^n]
n_reg               = 0.5
L_port              = 0.4           # [m]
A_port              = 0.0019634954  # [m^2]

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
