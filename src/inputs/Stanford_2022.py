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
sim_time = 12           # [s]
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
m_ox        = 10.0              # [kg]
P_ox_tank   = 4067908.4         # [Pa]
V_ox_tank   = 0.01289902263     # [m^3]
diam_out    = 0.127             # [m]
diam_in     = 0.12065           # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]
volume_err_tol = 0.01
P_dot_err_tol  = None

# ------------------------
# Ox Injector parameters
# ------------------------2
Cd_inj_ox   = 0.6#1.0
A_inj_ox = 27.08e-6      # [m^2]

# ------------------------
# Fuel tank & pressurant
# ------------------------
m_fuel       = 1.25
m_pres       = 0.092
P_pres_tank  = 4067908.4    # [Pa]
V_fuel_tank  = 0.00352530695  # [m^3]
V_pres_tank  = 0             # [m^3]  #NOTE: NO EXT PRESSURANT TANK!
diam_out_fuel = 0.127           # [m]
diam_in_fuel  = 0.12065           # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]


# ------------------------
# Fuel Injector parameters
# ------------------------
Cd_inj_fuel   = 1.0
A_inj_fuel = 1.2e-6 # [m^2] #NOTE: 2.4e-6 reported but also shown in a graph it seems like it corresponds to a peak discharge and 1.2 is more 


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
validation_files={
        #"P_cc":         "src/inputs/Stanford_2022/Stanford_2022_CC_Pressure.csv", #NOTE: empty csv
        "P_ox_tank":    "src/inputs/Stanford_2022/Stanford_2022_Ox_Tank_Pressure.csv",
        "P_fuel_tank":  "src/inputs/Stanford_2022/Stanford_2022_Fuel_Tank_Pressure.csv",
        "thrust":       "src/inputs/Stanford_2022/Stanford_2022_Thrust.csv"
}

