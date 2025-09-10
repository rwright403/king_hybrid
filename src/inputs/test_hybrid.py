# ------------------------
# Std Engine Input File (Hybrid)
# ------------------------

# Models
ox_tank_model = 2   # 1 = equilibrium, 2 = non-equilibrium
ox_inj_model  = 1   # injector model selection
cc_model      = 1   # hybrid chamber (fuel regression)
nozzle_model  = 1   # nozzle model

# Global settings
thrust_curve_graphs = True
mode = "full_stack"   # "ox_tank", "fuel_tank", or "full_stack"
save_path = None

# Environment
timestep = 1e-2       # [s]
sim_time = 3.0        # [s]
P_atm    = 101325     # [Pa]
T_atm    = 298.15     # [K]
rho_atm  = 1.225      # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------
fuel_name       = "HTPB"
fuel_properties = "HTPB C=0.9 H=0.1"   # required for hybrid solid fuel
oxidizer_name   = "N2O"

# ------------------------
# Tank parameters (oxidizer tank only, no liquid fuel tank here)
# ------------------------
V_tank      = 0.010     # [m^3] oxidizer tank volume (10 L)
diam_out    = 0.152     # [m] outer diameter (6 in)
diam_in     = 0.146     # [m] inner diameter (~5.75 in)
rho_wall    = 2700.0    # [kg/m^3] Al 6061 density
k_w         = 167.0     # [W/m-K] Al 6061 conductivity

# Oxidizer tank initial state
m_ox        = 5.0       # [kg] nitrous oxide mass
P_ox_tank   = 5e6       # [Pa] initial oxidizer tank pressure (~50 bar)
P_cc        = 3e6       # [Pa] chamber back pressure (~30 bar)
volume_err_tol = 1e-6   # [-] convergence tolerance for volume solver
P_dot_err_tol  = 1e-3   # [-] convergence tolerance for dP/dt match

# ------------------------
# Injector parameters
# ------------------------
Cd_inj   = 0.85           # [-] discharge coefficient
A_inj_ox = 1e-6           # [m^2] oxidizer injector area

# ------------------------
# Chamber (hybrid regression model)
# ------------------------
L_star    = 20            # [m]
m_fuel_i  = 2.0           # [kg] initial solid fuel mass
rho_fuel  = 930.0         # [kg/m^3] solid fuel density
a_reg     = 0.2e-3        # regression rate coefficient
n_reg     = 0.5           # regression rate exponent
L_port    = 0.5           # [m] port length
A_port_i  = 1e-4          # [m^2] initial port cross-sectional area

# ------------------------
# Nozzle
# ------------------------
d_throat = 0.0381         # [m] throat diameter
expratio = 2.0            # A_exit / A_throat
