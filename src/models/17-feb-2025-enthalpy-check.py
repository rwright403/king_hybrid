import CoolProp.CoolProp as CP
from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical
import matplotlib.pyplot as plt
import traceback

# Fluid
fluid = 'N2O'

# Global Constants:
R_U = 8.31446 #J/(mol K)
T_REF = 298.15 #K

n2o_global = Chemical('nitrous oxide', T=T_REF)
PC = n2o_global.Pc
TC = n2o_global.Tc
OMEGA = n2o_global.omega

MW = (n2o_global.MW/1000) #n2o.MW in global/mol --> converted to kglobal/mol
KAPPA = 0.37464 + 1.5422*n2o_global.omega - 0.26992*n2o_global.omega**2
b = 0.07780*(R_U*TC/PC)
g = 9.81 #m/s^2

# To shift the reference point, subtract h_new_ref from future calculations
def convert_enthalpy_to_nist_convention(h, T, P):
    T_ref_nist = 184.68  # K
    P_ref_nist = 101325 # Pa
    
    n2o_ig_nist = Chemical('N2O', T=T_ref_nist) 
    preos_nist = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_ref_nist, P=P_ref_nist)
    h_dep_nist = preos_nist.H_dep_g/MW # J/kg
    h_ideal_nist = n2o_ig_nist.H
    h_ref_nist = h_ideal_nist + h_dep_nist
    print("h_ref_nist: ", h_ref_nist )

    return h - h_ref_nist

def convert_int_energy_to_nist_convention(u, T, P):
    T_ref_nist = 184.68  # K
    P_ref_nist = 101325 # Pa
    
    n2o_ig_nist = Chemical('N2O', T=T_ref_nist) 
    preos_nist = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_ref_nist, P=P_ref_nist)
    u_dep_nist = preos_nist.H_dep_g/MW # J/kg
    u_ideal_nist = n2o_ig_nist.H
    u_ref_nist = u_ideal_nist + u_dep_nist
    print("u_ref_nist: ", u_ref_nist )

    return u - u_ref_nist



T_1 = 290 #K
P_1 = 2.5e6 #Pa
h_NIST = 416.64e3 #J/kg

"""
T_1 = 184.68
P_1 = 101325
h_NIST = 0.0067


T_1 = 270
P_1 = 2e6
h_NIST = 417.21e3


T_1 = 290 #K
P_1 = 2e6 #Pa
h_NIST = 439.63e3 #J/kg

T_1 = 270 #K
P_1 = 2e6#Pa
h_NIST = 417.21e3 #J/kg

T_1 = 300 #K
P_1 = 5.5e6 #Pa
h_NIST = 388.88E3 #J/kg
"""


### solve enthalpy with 
n2o_ig = Chemical('N2O', T=T_1) 
preos = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_1, P=P_1)
h_dep = preos.H_dep_g/MW # J/kg
h_ideal = n2o_ig.H #J/kg
h_1 = h_ideal + h_dep

print("h_1 before converting to nist convention: ", h_1)
print("h_1 components! ", h_1, " = ", h_ideal, " + ", h_dep)

h_1 = convert_enthalpy_to_nist_convention(h_1, T_1, P_1)

print("checking enthalpies: ", h_1, h_NIST, "difference: ", h_1 - h_NIST)
print("note phase: ", preos.phase)


h_cool = CP.PropsSI('H', 'T', T_1, 'P', P_1, fluid)
print("coolprop: ", h_cool)


u_dep = preos.U_dep_g/MW # J/kg
u_ideal = n2o_ig.U #J/kg
u_1 = u_ideal + u_dep


u_cool = CP.PropsSI('U', 'T', T_1, 'P', P_1, fluid)

print("checking internal energy: ", u_1, u_cool)
u_1 = convert_int_energy_to_nist_convention(u_1, T_1, P_1)
print("converted u_1 ", u_1)
"""
n2o_ig_g = Chemical('N2O', T=T_1, P=P_1) 
preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_1, P=P_1)
h_gas = (preos_g.H_dep_g/MW +n2o_ig_g.Cpg*(T_1 - T_REF))  

print("how the script is solving: ", h_gas)
"""


