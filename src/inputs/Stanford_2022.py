# ------------------------
# Std Engine Input File (Converted Test Case - Liquid Rocket)
# ------------------------
import numpy as np

# Models
ox_tank_model   = 2  
fuel_tank_model = 1
ox_inj_model    = 3
cc_model        = 2 
nozzle_model    = 1

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 0.005       # [s]
sim_time = 6           # [s]
P_atm    = 101325     # [Pa]
T_atm    = 284.0      # [K]
rho_atm  = 1.225      # [kg/m^3]

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
P_ox_tank   = 4067908.4         # [Pa]
V_tank      = 0.01289902263     # [m^3]
diam_out    = 0.127             # [m]
diam_in     = 0.12065           # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]
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
V_tank       = 0.0021875 # [m^3]
V_pres_tank  = 0.00133780695             # [m^3] #NOTE: STANFORD DIRECTLY PRESSURIZES FUEL TANK, V_tank is psace fuel takes up and V_pres_tank is the ullage vol of the same tank (this is defined because of how my program takes in inputs)
diam_out_fuel = 0.127           # [m]
diam_in_fuel  = 0.12065           # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]

# ------------------------
# Chamber (liquid engine)
# ------------------------
V_cc = 0.00141961201   # [m^3]

# ------------------------
# Nozzle
# ------------------------
d_throat = 0.0254
expratio = 3

# ------------------------
# Experimental validation data
# ------------------------
exp_thrust_file_path    = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_Thrust.csv'
exp_p_cc_file_path      = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_CC_Pressure.csv'
exp_p_ox_tank_file_path = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_Ox_Tank_Pressure.csv'
exp_p_fuel_tank_file_path = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_Fuel_Tank_Pressure.csv'

