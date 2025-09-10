# ------------------------
# Std Engine Input File (Converted Test Case - Liquid Rocket)
# ------------------------

# Models
ox_tank_model   = 2   # bens_ox_tank (non-equilibrium style)
fuel_tank_model = 4   # piping_real_fuel_tank_data
ox_inj_model    = 4   # injector model (from inj_model=4)
cc_model        = 2   # adiabatic_lre_cc
nozzle_model    = 1   # nozzle model

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
oxidizer_name   = "N2O"
fuel_name       = "Ethanol"
fuel_properties = None      # not needed for LRE
fuel_str        = "Ethanol"
pres_str        = "N2"      # pressurant gas

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 12.0
P_ox_tank   = 5.51581e6     # [Pa]
P_cc        = None          # not specified
V_ox_tank   = 0.25*np.pi*((0.0254*5.625)**2 - (0.0254*2.84)**2) * (0.0254*64)
diam_out    = None
diam_in     = None
rho_wall    = None
k_w         = None
volume_err_tol = 0.01
P_dot_err_tol  = None

# ------------------------
# Injector parameters
# ------------------------
Cd_inj   = 0.8      # guessed
A_inj_ox = 0.25*np.pi*(1.2*(0.0254*0.06793)**2)*6
A_inj_fuel = 0.25*np.pi*((0.0254*0.06793)**2)*5*6

# ------------------------
# Fuel tank (piping_real_fuel_tank_data)
# ------------------------
m_fuel       = None
m_pres       = None
P_pres_tank  = None
V_tank       = None
V_pres_tank  = None
diam_out_fuel = None
diam_in_fuel  = 0.0254*2.84
rho_wall_fuel = None
k_w_fuel      = None
T_tank        = T_atm
fuel_tank_pressure_filepath = r'./src/inputs/liquid_validation_data/Waterloo_2024/Waterloo_2024_Fuel_Tank_Pressure.csv'

# ------------------------
# Chamber (liquid engine)
# ------------------------
L_star    = None

# ------------------------
# Nozzle
# ------------------------
A_throat = 0.25*np.pi*(0.0254*(1.975))**2
A_exit   = A_throat * 4.1
d_throat = None
expratio = None

# ------------------------
# Experimental validation data
# ------------------------
exp_thrust_file_path    = r'./src/inputs/liquid_validation_data/Waterloo_2024/Waterloo_2024_Thrust.csv'
exp_p_cc_file_path      = r'./src/inputs/liquid_validation_data/Waterloo_2024/Waterloo_2024_CC_Pressure.csv'
exp_p_ox_tank_file_path = r'./src/inputs/liquid_validation_data/Waterloo_2024/Waterloo_2024_Oxidizer_Tank_Pressure.csv'
exp_p_fuel_tank_file_path = r'./src/inputs/liquid_validation_data/Waterloo_2024/Waterloo_2024_Fuel_Tank_Pressure.csv'
