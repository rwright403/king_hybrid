# ------------------------
# Std Engine Input File (Converted Test Case - Liquid Rocket)
# ------------------------

# Models
ox_tank_model   = 2   # bens_ox_tank (non-equilibrium style)
fuel_tank_model = 3   # simpleAdiabaticPressurizedTank
ox_inj_model    = 4   # injector model (from inj_model=4)
cc_model        = 2   # adiabatic_lre_cc
nozzle_model    = 1   # nozzle model

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 0.05       # [s]
sim_time = 7.0        # [s]
P_atm    = 101325     # [Pa]
T_atm    = 293.15     # [K]
rho_atm  = None       # not provided

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
m_ox        = 5.73815484158   # [kg]
P_ox_tank   = 5.2e6           # [Pa]
P_cc        = None            # not specified
V_tank      = 0.01177057403   # [m^3]
diam_out    = None
diam_in     = None
rho_wall    = None
k_w         = None
volume_err_tol = 0.01
P_dot_err_tol  = None

# ------------------------
# Injector parameters
# ------------------------
Cd_inj   = 0.6
A_inj_ox = 2.5 * 3.090605599220321e-5   # [m^2]
A_inj_fuel = 6.566075013100621e-6       # [m^2]

# ------------------------
# Fuel tank & pressurant
# ------------------------
m_fuel       = 1.14763096832    # [kg]
m_pres       = 0.12             # [kg] pressurant mass (estimated)
P_pres_tank  = 5.2e6            # [Pa] fuel tank pressure
V_tank       = 0.0037961342     # [m^3] fuel tank volume
V_pres_tank  = None             # not provided
diam_out_fuel = 0.0254 * 5      # [m]
diam_in_fuel  = None
rho_wall_fuel = None
k_w_fuel      = None

# ------------------------
# Chamber (liquid engine)
# ------------------------
L_star    = None

# ------------------------
# Nozzle
# ------------------------
A_throat = 0.001342491     # [m^2]
A_exit   = 0.0055046156    # [m^2]
d_throat = None
expratio = None

# ------------------------
# Experimental validation data
# ------------------------
exp_thrust_file_path    = None
exp_p_cc_file_path      = None
exp_p_ox_tank_file_path = None
exp_p_fuel_tank_file_path = None

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
