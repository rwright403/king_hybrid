# ------------------------
# Tomacz Test Case Input File
# ------------------------

# Models
ox_tank_model = 2     # non-equilibrium tank
ox_inj_model  = 3     
cc_model      = 3     # pass in downstream chamber pressures
nozzle_model  = None  # no nozzle

# Global settings
thrust_curve_graphs = True
mode = "ox_tank"
save_path = None

# Environment
timestep = 0.05         # [s]
sim_time = 4.9            # [s]  
P_atm    = 101325       # [Pa]
T_atm    = 296.0        # [K]  (approx lab ambient, adjust if needed)
rho_atm  = 1.225        # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------
oxidizer_name = "N2O"
fuel_name     = None
fuel_properties = None

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 0.180          # [kg] from Tomacz run
P_ox_tank   = 5.2e6          # [Pa]
V_ox_tank      = 0.25 * 3.14159 * (0.040**2) * 0.220   # [m^3] from cylinder geometry
diam_in     = 0.040          # [m]
diam_out    = 0.0475         # [m]
rho_wall    = 2770           # [kg/m^3] Al6061
k_w         = 237            # [W/m-K]
volume_err_tol = 0.01        # default
P_dot_err_tol  = 1e-3        # small tolerance for secant solve

# ------------------------
# Injector parameters
# ------------------------
Cd_inj_ox   = 0.45#0.57#0.45
A_inj_ox = 1.76714E-06  # [m^2]
A_inj_fuel = None  # no fuel injector

# ------------------------
# Chamber
# ------------------------
L_star   = None
m_fuel_i = None
rho_fuel = None
a_reg    = None
n_reg    = None
L_port   = None
A_port   = None

# ------------------------
# Nozzle
# ------------------------
d_throat = None
expratio = None

# ------------------------
# Validation Files
# ------------------------
validation_files = {
    "P_cc": "src/inputs/tomasz_test_case/tomasz_test_case_CC_Pressure.csv",
    "P_ox_tank": "src/inputs/tomasz_test_case/tomasz_test_case_Ox_Tank_Pressure.csv",
    "m_ox_tank": "src/inputs/tomasz_test_case/tomasz_test_case_m.csv",
    "m_dot_ox_tank": "src/inputs/tomasz_test_case/tomasz_test_case_m_dot.csv"
}
