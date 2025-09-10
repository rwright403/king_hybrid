# ------------------------
# Std Engine Input File (Converted Test Case)
# ------------------------

# Models
ox_tank_model = 1   # bens_ox_tank
ox_inj_model  = 1   # SPI injector (from inj_model=1)
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
T_atm    = 295.0      # [K]
rho_atm  = None       # not provided

# ------------------------
# Propellant properties
# ------------------------
fuel_name       = "paraffin"
fuel_properties = """#fuel paraffin  C 20   H 42    wt%=100.00
#h,KJ/mol=-13313.49  t(k)=298.15   rho,kg/m3=900
"""
oxidizer_name   = "N2O"

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 8.2           # [kg]
P_ox_tank   = 3920207.656   # [Pa]
P_cc        = None          # not specified
V_tank      = 0.012         # [m^3]
diam_out    = None
diam_in     = None
rho_wall    = None
k_w         = None
volume_err_tol = 0.01       # from all_error
P_dot_err_tol  = None

# ------------------------
# Injector parameters
# ------------------------
Cd_inj   = 1.0
A_inj_ox = 2.94229787e-5   # [m^2]
A_inj_fuel = None          # hybrid has no fuel injector

# ------------------------
# Chamber (hybrid regression model)
# ------------------------
L_star    = None
m_fuel_i  = 1.5          # [kg]
rho_fuel  = 900.0        # [kg/m^3]
a_reg     = 0.155e-3     # [m/s*(kg/s)^n]
n_reg     = 0.5
L_port    = 0.3852333    # [m]
A_port_i  = 0.0038319753 # [m^2]

# ------------------------
# Nozzle
# ------------------------
d_throat = None          # not directly provided
expratio = None          # not directly provided

# ------------------------
# Experimental validation data (kept for reference)
# ------------------------
exp_thrust_file_path   = r'./src/inputs/bens_validation_data/UofT_Deliverance_II/UofT_Deliverance_II_Thrust.csv'
exp_p_cc_file_path     = r'./src/inputs/bens_validation_data/UofT_Deliverance_II/UofT_Deliverance_II_CC_Pressure.csv'
exp_p_ox_tank_file_path = r'./src/inputs/bens_validation_data/UofT_Deliverance_II/UofT_Deliverance_II_Tank_Pressure.csv'
# exp_p_fuel_tank_file_path = None
