# ------------------------
# Std Engine Input File (Converted Test Case - Liquid Rocket)
# ------------------------
import numpy as np

# Models
ox_tank_model   = 2   # bens_ox_tank (non-equilibrium style)
fuel_tank_model = 1
ox_inj_model    = 3
cc_model        = 2   # adiabatic_lre_cc
nozzle_model    = 1   # nozzle model

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 0.05       # [s]
sim_time = 15.0       # [s]
P_atm    = 101325     # [Pa]
T_atm    = 293.15     # [K]
rho_atm  = None       # not provided

# ------------------------
# Propellant properties
# ------------------------
oxidizer_name   = "N2O"
fuel_name       = "Kerosene"      # chamber fuel
fuel_properties = None            # not needed for LRE
fuel_str        = "n-Dodecane"    # tank fuel string for CoolProp
pres_str        = "N2"            # pressurant gas

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 10.0
P_ox_tank   = 6894.76 * 590   # [Pa]
P_cc        = None            # not provided
V_tank      = (0.0254*42.47) * 0.25*np.pi*((0.0254*4.75)**2)  # [m^3]
diam_out    = None
diam_in     = None
rho_wall    = None
k_w         = None
volume_err_tol = 0.01
P_dot_err_tol  = None

# ------------------------
# Injector parameters
# ------------------------
Cd_inj   = 1.0
A_inj_ox = 27.08e-6      # [m^2]
A_inj_fuel = 2.4e-6      # [m^2]

# ------------------------
# Fuel tank & pressurant
# ------------------------
m_fuel       = 1.75
m_pres       = 0.12
P_pres_tank  = 6894.76 * 610   # [Pa]
V_tank       = (0.0254*13.90) * 0.25*np.pi*((0.0254*4.75)**2)  # [m^3]
V_pres_tank  = None
diam_out_fuel = 0.0254*4.75
diam_in_fuel  = None
rho_wall_fuel = None
k_w_fuel      = None

# ------------------------
# Chamber (liquid engine)
# ------------------------
L_star    = None  # not provided

# ------------------------
# Nozzle
# ------------------------
# Provided A_throat, A_exit but leave as None unless we compute
A_throat = 0.25*np.pi*((0.0254*1)**2)
A_exit   = 3*A_throat
d_throat = None
expratio = None

# ------------------------
# Experimental validation data
# ------------------------
exp_thrust_file_path    = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_Thrust.csv'
exp_p_cc_file_path      = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_CC_Pressure.csv'
exp_p_ox_tank_file_path = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_Ox_Tank_Pressure.csv'
exp_p_fuel_tank_file_path = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_Fuel_Tank_Pressure.csv'

# ------------------------
# Sensitivity Analysis (not part of runtime input, kept for reference)
# ------------------------
test_var_name = "Cd_inj"
min_bound = 0.5
max_bound = 0.7
num_iterations = 3

# ------------------------
# Sizing wizard estimates (kept for reference, not used in runtime)
# ------------------------
min_TW_ratio = 11
Cd_est = 0.6
mass_fraction_estimate = 0.2657741984
characteristic_len = 1.0
