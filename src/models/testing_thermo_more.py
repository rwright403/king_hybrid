from thermo.eos import PR, IG
from thermo import Chemical
import CoolProp.CoolProp as CP

# Define inputs
T = 290  # Temperature in Kelvin
P = 5e6  # Pressure in Pascals (5 MPa)

# Critical properties for the fluid (example: nitrous oxide)
Tc = 309.57  # Critical temperature in Kelvin
Pc = 7.257e6  # Critical pressure in Pascals
omega = 0.04  # Acentric factor

molar_mass = 44.013/1000 # in kg/mol, for N2O it's ~44.013 g/mol

#check phase with coolprop:
#print(CP.PhaseSI('T', T, 'P', P, "N2O" ))
# Initialize Peng-Robinson EOS
eos = IG(T=T, P=P)#, Tc=Tc, Pc=Pc, omega=omega)
print(eos.phase)


#rho_vap = 1/eos.V_g * molar_mass  # Get the specific volume (m^3/mol)
n2o = Chemical('nitrous oxide')
rho_liq = 1/eos.V_l_sat(T)  * molar_mass # Get the specific volume (m^3/mol)
    


print(f"At T = {T} K and P = {P / 1e6} MPa:")
print(f"Liquid density: {rho_liq:.2f} kg/m³")
#print(f"Vapor density: {rho_vap:.2f} kg/m³")
