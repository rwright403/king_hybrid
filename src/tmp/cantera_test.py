import numpy as np
import cantera as ct
from rocketcea.cea_obj_w_units import CEA_Obj
from src.tmp.cantera_thermo_helper import *
import cantera as ct


### Start: use CEA to obtain relevant species for Cantera
C = CEA_Obj(
    oxName="N2O",
    fuelName="C2H5OH",   # ethanol
    pressure_units="Pa",
    isp_units="sec",
    cstar_units="m/s",
    temperature_units="K",
    sonic_velocity_units="m/s",
    enthalpy_units="J/kg",
    density_units="kg/m^3",
    specific_heat_units="J/kg-K",
)

P_cc = 3e6
OF   = 6.0
expratio = 40.0

molWtD, massFracD = C.get_SpeciesMassFractions(Pc=P_cc, MR=OF, eps=expratio)


print("Species in RocketCEA equilibrium:")
for k, v in massFracD.items():
    print(k, f"{v[1]*100:.2f}")

T_cc = C.get_Tcomb(Pc=P_cc, MR=OF)
MW, gamma = C.get_Chamber_MolWt_gamma(Pc=P_cc, MR=OF, eps=expratio)
cstar_cea = C.get_Cstar(Pc=P_cc, MR=OF)

print("MW: ", MW)
print("gamma: ", gamma)
print("cstar: ", cstar_cea)
print("T_cc: ", T_cc, "[K]")





"""
NOTE: in rocket cea massFracD values are ordered by station like:
[ injection face , chamber , throat , exit ]
                    ^
                    |
                    1st idx
"""
      


# load NASA database
full_species = {S.name: S for S in ct.Species.list_from_file("nasa_gas.yaml")}

# map RocketCEA species to Cantera names (RocketCEA uses * for atoms sometimes)
def map_species_name(name):
    return name.replace("*", "")  # '*O' -> 'O', '*N2' -> 'N2', etc.

rocketcea_species = [map_species_name(sp) for sp in massFracD.keys()]

# add ethanol + nitrous oxide explicitly
rocketcea_species += ["C2H5OH", "N2O"]

species = [full_species[s] for s in set(rocketcea_species) if s in full_species]

gas = ct.Solution(thermo="ideal-gas", kinetics="none", species=species)


# Convert O/F into mole fractions
MW_fuel = full_species["C2H5OH"].molecular_weight
MW_ox   = full_species["N2O"].molecular_weight

n_fuel = 1.0
m_fuel = n_fuel * MW_fuel
m_ox   = OF * m_fuel
n_ox   = m_ox / MW_ox

X = {"C2H5OH": n_fuel, "N2O": n_ox}

gas.TPX = 298.0, P_cc, X
gas.equilibrate("HP") #const enthalpy and pressure process

print(f"T [K]   = {gas.T:.2f}")
print(f"P [bar] = {gas.P/1e5:.2f}")
print("Equilibrium composition (mole %):")
for sp, frac in zip(gas.species_names, gas.X*100):
    if frac > 1e-2:
        print(f"{sp:>6s}: {frac:7.3f}")

### NOTE: seeing some species mass fractions off by 7% between cantera and CEA

gamma_ct = gas.cp_mole / gas.cv_mole
print(f'Cantera specific heat ratio: {gamma_ct: .4f}')




derivs = get_thermo_derivatives(gas)

dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(
    gas, derivs[0], derivs[1], derivs[2]
    )

print(f'Cp = {cp: .2f} J/(K kg)')

print(f'(d log V/d log P)_T = {dlogV_dlogP_T: .4f}')
print(f'(d log V/d log T)_P = {dlogV_dlogT_P: .4f}')

print(f'gamma_s = {gamma_s: .4f}')


print("mol wt, R_spec: ", gas.mean_molecular_weight, (ct.gas_constant/gas.mean_molecular_weight))
speed_sound = np.sqrt(ct.gas_constant * gas.T * gamma_s / gas.mean_molecular_weight)
print(f'Speed of sound = {speed_sound: .1f} m/s')

def calculate_c_star(gamma, temperature, molecular_weight):
    return (
        np.sqrt(ct.gas_constant * temperature / (molecular_weight * gamma)) *
        np.power(2 / (gamma + 1), -(gamma + 1) / (2*(gamma - 1)))
        )


entropy_chamber = gas.s
enthalpy_chamber = gas.enthalpy_mass
mole_fractions_chamber = gas.X
gamma_chamber = gamma_s

c_star = calculate_c_star(gamma_chamber, gas.T, gas.mean_molecular_weight)
print(f'c-star: {c_star: .1f} m/s')
print('Error in c-star: '
      f'{100*np.abs(c_star - cstar_cea)/cstar_cea: .3e} %'
      )