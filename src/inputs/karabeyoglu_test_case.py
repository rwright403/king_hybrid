# ------------------------
# Karabeyoglu Test Case Input File
# ------------------------

# Models
ox_tank_model = 2     # non-equilibrium tank
ox_inj_model  = 1     # SPI/orifice model
cc_model      = 3     # validation downstream Pc
nozzle_model  = None  # no nozzle

# Global settings
thrust_curve_graphs = True
mode = "ox_tank"
save_path = None

# Environment
#TODO: PAPER REPORTED THEIR TIMESTEP WAS 1E-3 BUT I THINK WE CAN GO HIGHER
timestep = 1e-3        # [s]
sim_time = 7.0         # [s]
P_atm    = 1.0e5       # [Pa]
T_atm    = 286.5       # [K]
rho_atm  = 1.225       # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------
oxidizer_name = "N2O"
fuel_name     = None
fuel_properties = None

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 20.0          # [kg]
P_ox_tank   = 4.5e6         # [Pa] initial
V_tank      = 0.0354        # [m^3]
diam_in     = 0.1905        # [m]
diam_out    = 0.1905 + (3.18e-3)*2  # [m]
rho_wall    = 2770          # [kg/m^3]
k_w         = 237           # [W/m-K]
volume_err_tol = 0.01
P_dot_err_tol  = 1e-3

# ------------------------
# Injector parameters
# ------------------------
Cd_inj   = 0.425
A_inj_ox = 1.0e-4    # [m^2] (guess from paper)
A_inj_fuel = None

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
    "P_ox_tank": "src/inputs/karabeyoglu_test_case/karabeyoglu_P_tank.csv",
    "T_liq_ox_tank": "src/inputs/karabeyoglu_test_case/karabeyoglu_T_liq.csv",
    "T_gas_ox_tank": "src/inputs/karabeyoglu_test_case/karabeyoglu_T_gas.csv",
    "T_sat_ox_tank": "src/inputs/karabeyoglu_test_case/karabeyoglu_T_sat.csv",
}
