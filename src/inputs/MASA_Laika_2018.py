# ------------------------
# Std Engine Input File (Converted Test Case - Liquid Rocket)
# ------------------------

# Models
ox_tank_model   = 2   # bens_ox_tank (non-equilibrium style)
fuel_tank_model = 3   # simpleAdiabaticPressurizedTank
ox_inj_model    = 1   # SPI injector (from inj_model=1)
cc_model        = 2   # adiabatic_lre_cc
nozzle_model    = 1   # nozzle model

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 0.05       # [s]
sim_time = 3.4        # [s]
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
m_ox        = 4.48         # [kg]
P_ox_tank   = 5.171e6      # [Pa]
P_cc        = None         # not given
V_tank      = 6.4e-3       # [m^3]
diam_out    = None
diam_in     = None
rho_wall    = None
k_w         = None
volume_err_tol = 0.01
P_dot_err_tol  = None

# ------------------------
# Injector parameters
# ------------------------
Cd_inj   = 0.66
A_inj_ox = 7.471705e-5     # [m^2]
A_inj_fuel = 1.36284e-5    # [m^2]

# ------------------------
# Fuel tank & pressurant
# ------------------------
m_fuel       = 1.12        # [kg]
m_pres       = 0.12        # [kg] pressurant mass (estimated)
P_pres_tank  = 4.82633e6   # [Pa] (from P_fueltank)
V_tank       = 2.16e-3     # [m^3] fuel tank volume
V_pres_tank  = None        # not given
diam_out_fuel = 0.0254*5   # [m] ID_PROPTANK, assuming as inner diameter
diam_in_fuel  = None
rho_wall_fuel = None
k_w_fuel      = None

# ------------------------
# Chamber (liquid engine)
# ------------------------
L_star    = None           # not given, can be filled if available

# ------------------------
# Nozzle
# ------------------------
# A_throat and A_exit provided → compute equivalent geometry
d_throat = None             # can compute if needed
expratio = None             # can compute if needed

