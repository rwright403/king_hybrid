# ------------------------
# Liquid Rocket Input File
# ------------------------

# Models
ox_tank_model   = 1
fuel_tank_model = 1
ox_inj_model    = 1
cc_model        = 2
nozzle_model    = 1

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
TIMESTEP = 1e-3
sim_time = 0.5
P_atm = 101325

# CEA object
from rocketcea.cea_obj_w_units import CEA_Obj

oxidizer_name = "N2O"
fuel_name     = "Ethanol"

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
    P_atm=P_atm,
    all_error=1e-3,
)

# Fuel Tank
from src.utils.enum import FillType  # wherever you defined FillType

fuel_tank_kwargs = dict(
    m_fuel=3.0,
    m_pres=0.5,                # [kg] pressurant mass
    P_pres_tank=2e7,           # [Pa] initial pressurant bottle pressure
    P_atm=P_atm,
    T_atm=298.0,               # [K] ambient
    rho_atm=1.225,             # [kg/m^3] ambient air density
    V_fuel_tank=0.01,         # [m^3]
    V_pres_tank=0.002,         # [m^3] pressurant bottle volume
    diam_out=0.12,             # [m] outer tank diameter
    diam_in=0.10,              # [m] inner tank diameter
    rho_wall=7800,             # [kg/m^3] (steel)
    k_w=15.0,                  # [W/m-K] wall conductivity
    Cd=0.85,
    A_inj=1e-5,                # [m^2] injector area
    fuel_str="Ethanol",
    pres_str="Nitrogen",       # or "Helium", depending on your pressurant
    fill_type=FillType.ISOTHERMAL,  # choose FillType.ISOTHERMAL if you want
)


# Oxidizer Injector
ox_inj_kwargs = dict(
    Cd=0.85,
    A_inj=5e-5,
)

# Fuel Injector
fuel_inj_kwargs = dict(
    Cd=0.85,
    A_inj=1e-5,
)

# Combustion Chamber
cc_kwargs = dict(
    L_star=20,
    C=C
)

# Nozzle
nozzle_kwargs = dict(
    d_throat=0.0381,
    expratio=2.0,
    P_atm=P_atm,
    C=C
)
