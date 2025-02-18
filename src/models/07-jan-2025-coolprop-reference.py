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
P_REF = 101325 #Pa

n2o_global = Chemical('nitrous oxide', T=T_REF)
PC = n2o_global.Pc
TC = n2o_global.Tc
OMEGA = n2o_global.omega

MW = (n2o_global.MW/1000) #n2o.MW in global/mol --> converted to kglobal/mol
KAPPA = 0.37464 + 1.5422*n2o_global.omega - 0.26992*n2o_global.omega**2
b = 0.07780*(R_U*TC/PC)
g = 9.81 #m/s^2

# Triple point (CoolProp reference)
T_triple = 182.33  # K
P_triple = 87910  # Pa

# To shift the reference point, subtract h_new_ref from future calculations
def convert_enthalpy_to_nist_convention(T, P):
    T_ref = T_triple # K
    P_ref = P_triple # Pa
    h_ref = CP.PropsSI('H', 'T', T_ref, 'P', P_ref, fluid)

    h = CP.PropsSI('H', 'T', T, 'P', P, fluid)
    return h - h_ref

# Example: Calculate enthalpy at T = 250 K and P = 1 atm with the adjusted reference
T_example = 290  # K
P_example = 2.5 # Pa
h_adjusted = convert_enthalpy_to_nist_convention(T_example, P_example)
print(f"Adjusted enthalpy at T = {T_example} K and P = {P_example} Pa: {h_adjusted} J/kg")
### checking thermo departure reference enthalpy!!!!


"""
#1 try critical point
h_1 = None
T_1 = TC
P_1 = PC
preos_1 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_1, P=P_1)
try:
    h_1 = preos_1.H_dep_g
except Exception as e:
    try:
        print("case 1 is liquid")
        h_1 = preos_1.H_dep_l
    except Exception as e:
        print("case 1 does not work?")
    

#2 try normal boiling point
h_2 = None
T_2 = 184.68
P_2 = 101325
preos_2 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_2, P=P_2)
try:
    h_2 = preos_2.H_dep_g
except Exception as e:
    try:
        print("case 2 is liquid")
        h_2 = preos_2.H_dep_g
    except Exception as e:
        print("case 2 does not work?")


#3 try atm same as ^
h_3 = None
T_3 = 298.15
P_3 = 101325
preos_3 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_3, P=P_3)
try:
    h_3 = preos_3.H_dep_g
except Exception as e:
    try:
        print("case 3 is liquid")
        h_3 = preos_3.H_dep_l
    except Exception as e:
        print("case 3 does not work?")


#4 try T at 0 degrees celsius?
h_4 = None
T_4 = 273.15
P_4 = 101325
preos_4 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_4, P=P_4)
try:
    h_4 = preos_4.H_dep_g
except Exception as e:
    try:
        print("case 4 is liquid")
        h_4 = preos_4.H_dep_l
    except Exception as e:
        print("case 4 does not work?")


#5try T at triple point
h_5 = None
T_5 = 184.68
P_5 = 87850
preos_5 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_5, P=P_5)
try:
    h_5 = preos_5.H_dep_g
except Exception as e:
    try:
        print("case 5 is liquid")
        h_5 = preos_5.H_dep_l
    except Exception as e:
        print("case 5 does not work?")

print("looking for enthalpy reference state: ", h_1, h_2, h_3, h_4, h_5)
"""