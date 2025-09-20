# ------------------------
# Results from the U of Washington Boundless IREC Final Report (second hotfire)
# ------------------------

# Models
ox_tank_model = 2   
ox_inj_model  = 1   
cc_model      = 1
nozzle_model  = 1

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 0.005       # [s]
sim_time = 18       # [s]
P_atm    = 101325     # [Pa]
T_atm    = 292.0      # [K]
rho_atm  = 1.225      # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------
fuel_name       = "paraffin"
fuel_properties = """#fuel paraffin  C 20   H 42    wt%=100.00
#h,KJ/mol=-13313.49  t(k)=298.15   rho,kg/m3=900
"""
#BUG: above commented out! check this!
oxidizer_name   = "N2O"

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 15.4221       # [kg]
P_ox_tank   = 4933902       # [Pa]
V_tank      = 0.02532569531   # [m^3]
diam_out    = 0.2032        # [m]
diam_in     = 0.1905        # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237
volume_err_tol = 0.01       # mapped from all_error
P_dot_err_tol  = None

# ------------------------
# Injector parameters
# ------------------------
Cd_inj   = 0.6
A_inj_ox = 8.75591e-5   # [m^2]
A_inj_fuel = None           # hybrid has no fuel injector

# ------------------------
# Chamber (hybrid regression model)
# ------------------------
L_star    = 38.08382787  # [m] solved from volume of cc (not including fuel grain) and throat area
m_fuel_i  = 5.44311      # [kg]
rho_fuel  = 900.0        # [kg/m^3]
a_reg     = 0.155e-3     # [m/s*(kg/s)^n]
n_reg     = 0.5
L_port    = 0.41130982   # [m]
A_port  = 0.001302041  # [m^2]

# ------------------------
# Nozzle
# ------------------------
d_throat = 0.04071621409    
expratio = 4.721060934

# ------------------------
# Validation Files
# ------------------------
validation_files={
        "P_cc": "src/inputs/UofW_Boundless/UofW_Boundless_CC_Pressure.csv",
        "P_ox_tank": "src/inputs/UofW_Boundless/UofW_Boundless_Tank_Pressure.csv",
        "thrust": "src/inputs/UofW_Boundless/UofW_Boundless_Thrust.csv"
    }

