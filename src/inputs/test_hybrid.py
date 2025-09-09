# ------------------------
# Hybrid Rocket Input File
# ------------------------

# Models
ox_tank_model = 2
ox_inj_model  = 1  
cc_model      = 1   
nozzle_model  = 1

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
TIMESTEP = 1e-3
sim_time = 2.0
P_atm = 101325


# CEA object
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel

fuel_name = "HTPB"
fuel_properties = "HTPB C=0.9 H=0.1"
oxidizer_name = "N2O"

# Register new solid fuel with RocketCEA
add_new_fuel(fuel_name, fuel_properties)

# Global CEA object
C = CEA_Obj(
    oxName=oxidizer_name,
    fuelName=fuel_name,
    pressure_units="Pa",
    isp_units="sec",
    cstar_units="m/s",
    temperature_units="K",
    sonic_velocity_units="m/s",
    enthalpy_units="kJ/kg",
    density_units="kg/m^3",
    specific_heat_units="kJ/kg-K",
)


# Equilibrium Oxidizer Tank
"""
ox_tank_kwargs = dict(
    oxidizer="N2O",
    m_ox=5.0,
    V_tank=0.01,
    P_tank=5e6,
    P_atm=101325,
    all_error=1e-3,
)
"""

# Equilibrium Oxidizer Tank
"""
ox_tank_kwargs = dict(
    oxidizer="N2O",
    m_ox=5.0,
    V_tank=0.01,
    P_tank=5e6,
    P_atm=101325,
    all_error=1e-3,
)
"""

# Non-Equilibrium Oxidizer Tank

ox_tank_kwargs = dict(
    m_nos=5.0,              # [kg] initial nitrous oxide mass
    P_tank=5e6,             # [Pa] initial tank pressure (~50 bar)
    P_cc=3e6,               # [Pa] chamber back pressure (~30 bar)
    P_atm=101325,           # [Pa] ambient pressure
    T_atm=298.15,           # [K] ambient temperature (~25 °C)
    rho_atm=1.225,          # [kg/m3] air density at STP
    V_tank=0.010,           # [m3] total tank volume (10 L)
    diam_out=0.152,         # [m] outer tank diameter (6 in)
    diam_in=0.146,          # [m] inner tank diameter (~5.75 in)
    rho_wall=2700.0,        # [kg/m3] Al 6061 density
    k_w=167.0,              # [W/m-K] Al 6061 thermal conductivity
    volume_err_tol=1e-6,    # [-] convergence tolerance for volume solver
    P_dot_err_tol=1e-3,     # [-] convergence tolerance for dP/dt match
)



# Oxidizer Injector
ox_inj_kwargs = dict(
    Cd=0.85,
    A_inj=1e-5,
)

# Hybrid Chamber (fuel regression)
cc_kwargs = dict(
    L_star = 20,
    m_fuel_i=2.0,
    rho_fuel=930.0,
    a=0.2e-3,
    n=0.5,
    L=0.5,
    A_port_i=1e-4,
    C=C
)

# Nozzle
nozzle_kwargs = dict(
    d_throat=0.0381,   # throat diameter [m]
    expratio=2.0,    # A_exit / A_throat
    P_atm=P_atm,
    C=C
)
