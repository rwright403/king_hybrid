# ------------------------
# Std Engine Input File
# ------------------------



# Models
ox_tank_model = 2   # 1 = equilibrium, 2 = non-equilibrium
ox_inj_model  = 1   # injector model selection
cc_model      = 1   # chamber model
nozzle_model  = 1   # nozzle model


# Global settings
thrust_curve_graphs = True

# ------------------------
# Simulation control
# ------------------------
# Options:
#   "ox_tank"   → simulate oxidizer tank only
#   "fuel_tank" → simulate fuel/pressurant tank only
#   "full_stack" → simulate full engine (default)
mode = "full_stack"
save_path = None

# Environment
TIMESTEP = 1e-2       # [s]
sim_time = 3.0        # [s]
P_atm    = 101325     # [Pa]
T_atm    = 298.15     # [K]
rho_atm  = 1.225      # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------
fuel_name       = "HTPB"
fuel_properties = "HTPB C=0.9 H=0.1"
oxidizer_name   = "N2O"
fuel_str        = "Ethanol"   # for liquid-fuel config
pres_str        = "Nitrogen"  # pressurant gas (alt: "Helium")

# ------------------------
# Tank parameters (both ox & fuel)
# ------------------------
V_tank      = 0.010     # [m^3] tank volume
V_pres_tank = 0.002     # [m^3] pressurant bottle volume
diam_out    = 0.152     # [m] outer diameter
diam_in     = 0.146     # [m] inner diameter
rho_wall    = 2700.0    # [kg/m^3] wall material density
k_w         = 167.0     # [W/m-K] wall conductivity

# Oxidizer tank initial state
m_nos       = 5.0       # [kg] nitrous oxide mass
P_ox_tank   = 5e6       # [Pa] initial oxidizer tank pressure
P_cc        = 3e6       # [Pa] chamber back pressure
volume_err_tol = 1e-6
P_dot_err_tol  = 1e-3

# Fuel tank & pressurant
m_fuel       = 3.0      # [kg]
m_pres       = 0.5      # [kg] pressurant mass
P_pres_tank  = 2e7      # [Pa] pressurant bottle pressure

# ------------------------
# Injector parameters
# ------------------------
Cd_inj = 0.85           # [-]
A_inj_ox = 1e-6         # [m^2] oxidizer injector area
A_inj_fuel = 1e-5       # [m^2] fuel injector area

# ------------------------
# Chamber (fuel regression)
# ------------------------
L_star     = 20         # [m]
m_fuel_i   = 2.0        # [kg] initial solid fuel
rho_fuel   = 930.0      # [kg/m^3]
a_reg      = 0.2e-3     # [m/s*(kg/s)^n] regression rate coeff
n_reg      = 0.5        # [-] regression rate exponent
L_port     = 0.5        # [m] port length
A_port_i   = 1e-4       # [m^2] initial port cross-sectional area

# ------------------------
# Nozzle
# ------------------------
d_throat   = 0.0381     # [m]
expratio   = 2.0        # A_exit / A_throat
