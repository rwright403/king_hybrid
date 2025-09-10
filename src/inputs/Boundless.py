# ------------------------
# Std Engine Input File (Converted Test Case)
# ------------------------

# Models
ox_tank_model = 1   # bens_ox_tank
ox_inj_model  = 1   # generic injector
cc_model      = 1   # hybrid_cc_w_fuel_grain
nozzle_model  = 1   # nozzle model

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 0.05       # [s]
sim_time = 14.5       # [s]
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
m_ox        = 14.61751413   # [kg]
P_ox_tank   = 4933902       # [Pa]
P_cc        = None          # not specified
V_tank      = 0.024956005   # [m^3]
diam_out    = None
diam_in     = None
rho_wall    = None
k_w         = None
volume_err_tol = 0.01       # mapped from all_error
P_dot_err_tol  = None

# ------------------------
# Injector parameters
# ------------------------
Cd_inj   = 1.0
A_inj_ox = 2.364094414e-5   # [m^2]
A_inj_fuel = None           # hybrid has no fuel injector

# ------------------------
# Chamber (hybrid regression model)
# ------------------------
L_star    = None
m_fuel_i  = 5.44311      # [kg]
rho_fuel  = 900.0        # [kg/m^3]
a_reg     = 0.155e-3     # [m/s*(kg/s)^n]
n_reg     = 0.5
L_port    = 0.41130982   # [m]
A_port_i  = 0.001302041  # [m^2]

# ------------------------
# Nozzle
# ------------------------
d_throat = None           # not directly provided
expratio = None           # not directly provided
