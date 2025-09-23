# ------------------------
# Std Engine Input File (Converted Test Case)
# ------------------------

# Models
ox_tank_model = 2   
ox_inj_model  = 3   
cc_model      = 1   # hybrid_cc_w_fuel_grain
nozzle_model  = 1   # nozzle model

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 0.005       # [s]
sim_time = 6       # [s]
P_atm    = 101325     # [Pa]
T_atm    = 284.0      # [K]
rho_atm  = 1.225      # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------
fuel_name       = "paraffin" # NOTE: had issues with the rocketcea add fuel method so I modified library and added this card
oxidizer_name   = "N2O"

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 8.2           # [kg]
P_ox_tank   = 3920207.656   # [Pa]
V_tank      = 0.012         # [m^3]
diam_out    = 0.1397        # [m]
diam_in     = 0.128524      # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]
volume_err_tol = 0.001       # from all_error
P_dot_err_tol = None #TODO: RM

# ------------------------
# Injector parameters
# ------------------------
Cd_inj   = 0.611 # UTAT theoretically determined cd published in their report
A_inj_ox = 6.54E-05   # [m^2]
A_inj_fuel = None          # hybrid has no fuel injector

# ------------------------
# Chamber (hybrid regression model)
# ------------------------
V_pre_post_cc       = 0.00118174778  # [m^2]
m_fuel_i            = 1.5          # [kg]
rho_fuel            = 900.0        # [kg/m^3]
a_reg               = 0.155e-3     # [m/s*(kg/s)^n]
n_reg               = 0.5
L_port              = 0.3852333    # [m]
A_port              = 0.0038319753 # [m^2]

# ------------------------
# Nozzle
# ------------------------
d_throat = 0.037084      # [m]
expratio = 4.894398574   # [m]

# ------------------------
# Validation Files
# ------------------------
validation_files={
        "P_cc": "src/inputs/UofT_Deliverance_II/UofT_Deliverance_II_CC_Pressure.csv",
        "P_ox_tank": "src/inputs/UofT_Deliverance_II/UofT_Deliverance_II_Tank_Pressure.csv",
        "thrust": "src/inputs/UofT_Deliverance_II/UofT_Deliverance_II_Thrust.csv"
}