# ------------------------
# MASA Laika 3/17/18 Hotfire Analysis
# ------------------------

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
sim_time = 5.0           # [s]
P_atm    = 101325     # [Pa]
T_atm    = 292.0      # [K]
rho_atm  = 1.225      # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------
oxidizer_name   = "N2O"
fuel_name       = "Ethanol"
fuel_properties = None     # not needed for LRE
fuel_str        = "Ethanol"
pres_str        = "N2"     # pressurant gas

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 4.48         # [kg]
P_ox_tank   = 5.171e6      # [Pa]
V_ox_tank   = 6.4e-3       # [m^3]
diam_out    = 0.135128      # [m]
diam_in     = 0.127         # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]
volume_err_tol = 0.01
P_dot_err_tol  = None

# ------------------------
# Ox Injector parameters
# ------------------------
Cd_inj_ox   = 0.6
A_inj_ox = 7.471705e-5     # [m^2]


# ------------------------
# Fuel tank & pressurant
# ------------------------
m_fuel       = 1.12        # [kg]
m_pres       = 0.092        # [kg] pressurant mass (estimated)
P_pres_tank  = 4.82633e6   # [Pa] (from P_fueltank)
V_fuel_tank  = 3.25E-03     # [m^3] fuel tank volume
V_pres_tank  = 0.0        # [m^3]  #NOTE: NO EXT PRESSURANT TANK!
diam_out_fuel = 0.135128      # [m]
diam_in_fuel  = 0.127       # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]


# ------------------------
# Fuel Injector parameters
# ------------------------
Cd_inj_fuel   = 0.6 #NOTE: GUESS
A_inj_fuel = 1.36284e-5    # [m^2]

# ------------------------
# Chamber (liquid engine)
# ------------------------
V_cc      = 0.0032221466077      # [m^3]

# ------------------------
# Nozzle
# ------------------------
d_throat = 0.0360426        # [m]
expratio = 5.168493033     

# ------------------------
# Experimental validation data
# ------------------------
validation_files={
        "P_cc": "src/inputs/MASA_Laika_2018/MASA_Laika_2018_CC_Pressure.csv",
        "P_ox_tank": "src/inputs/MASA_Laika_2018/MASA_Laika_2018_Ox_Tank_Pressure.csv",
        "P_fuel_tank": "src/inputs/MASA_Laika_2018/MASA_Laika_2018_Fuel_Tank_Pressure.csv",
        "thrust": "src/inputs/MASA_Laika_2018/MASA_Laika_2018_Thrust.csv"
}