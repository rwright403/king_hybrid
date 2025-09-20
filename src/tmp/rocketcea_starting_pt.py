from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel

# ----------------------------------------
# 1. Define paraffin fuel (C32H66 surrogate)
# ----------------------------------------
# Format: "name   elements ... wt%=100"
paraffin_card = """paraffin  C 32  H 66   wt%=100.0"""

# Register new fuel
add_new_fuel("PARAFFIN", paraffin_card) #NOTE: ASSUMING C 32 H 66



C = CEA_Obj(
    oxName="N2O",
    fuelName="PARAFFIN",
    pressure_units="Pa",
    isp_units="sec",
    cstar_units="m/s",
    temperature_units="K",
    sonic_velocity_units="m/s",
    enthalpy_units="kJ/kg",
    density_units="kg/m^3",
    specific_heat_units="kJ/kg-K",
)

# ----------------------------------------
# 3. Define conditions
# ----------------------------------------
P_cc = 3e6         # chamber pressure [bar]
OF = 6.0               # O/F ratio
expratio = 40.0       # nozzle expansion ratio

# ----------------------------------------
# 4. Query RocketCEA for equilibrium properties
# ----------------------------------------
species_mass_fract = C.get_SpeciesMassFractions(Pc=P_cc, MR=OF, eps=expratio)
for i in species_mass_fract:
    print(i, "\n")

