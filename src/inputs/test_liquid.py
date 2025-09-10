# ------------------------
# Std Engine Input File
# ------------------------

# Models
ox_tank_model   = 2   # 1 = equilibrium, 2 = non-equilibrium
fuel_tank_model = 1   # pressurant-driven liquid tank
ox_inj_model    = 1   # injector model selection
cc_model        = 2   # chamber model (liquid engine)
nozzle_model    = 1   # nozzle model

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 1e-3       # [s]
sim_time = 0.5        # [s]
P_atm    = 101325     # [Pa]
T_atm    = 298.15     # [K] (~25 °C)
rho_atm  = 1.225      # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------
oxidizer_name   = "N2O"
fuel_name       = "Ethanol"
fuel_properties = None         # NOTE: only for defining custom fuel ie hybrid, set to None for lre
fuel_str        = "Ethanol"    # for fuel tank kwargs
pres_str        = "Nitrogen"   # pressurant gas (alt: "Helium")

# ------------------------
# Tank parameters (both ox & fuel)
# ------------------------
V_tank      = 0.010     # [m^3] oxidizer tank volume
V_pres_tank = 0.002     # [m^3] pressurant bottle volume
diam_out    = 0.152     # [m] outer diameter (ox tank default)
diam_in     = 0.146     # [m] inner diameter (ox tank default)
rho_wall    = 2700.0    # [kg/m^3] Al 6061 density
k_w         = 167.0     # [W/m-K] Al 6061 conductivity

# Oxidizer tank initial state
m_ox       = 5.0       # [kg] nitrous oxide mass
P_ox_tank   = 5e6       # [Pa] initial oxidizer tank pressure
P_cc        = 3e6       # [Pa] chamber back pressure
volume_err_tol = 1e-6
P_dot_err_tol  = 1e-3

# Fuel tank & pressurant
m_fuel       = 3.0      # [kg]
m_pres       = 0.5      # [kg] pressurant mass
P_pres_tank  = 2e7      # [Pa] pressurant bottle pressure
# geometry differs from ox tank:
diam_out_fuel = 0.12    # [m] outer diameter
diam_in_fuel  = 0.10    # [m] inner diameter
rho_wall_fuel = 7800    # [kg/m^3] steel
k_w_fuel      = 15.0    # [W/m-K] steel conductivity

# ------------------------
# Injector parameters
# ------------------------
Cd_inj    = 0.85           # [-]
A_inj_ox  = 5e-5           # [m^2] oxidizer injector area
A_inj_fuel = 1e-5          # [m^2] fuel injector area

# ------------------------
# Chamber (liquid engine)
# ------------------------
L_star   = 20              # [m]
# (liquid chamber → no regression model, so no m_fuel_i, a_reg, n_reg, etc.)

# ------------------------
# Nozzle
# ------------------------
d_throat = 0.0381          # [m]
expratio = 2.0             # A_exit / A_throat
