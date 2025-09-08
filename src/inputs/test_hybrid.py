# ------------------------
# Hybrid Rocket Input File
# ------------------------

# Models
ox_tank_model = 1
inj_model     = 1  
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


# Oxidizer Tank
ox_tank_kwargs = dict(
    oxidizer="N2O",
    m_ox=5.0,
    V_tank=0.01,
    P_tank=5e6,
    P_atm=101325,
    all_error=1e-3,
)

# Oxidizer Injector
inj_kwargs = dict(
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
