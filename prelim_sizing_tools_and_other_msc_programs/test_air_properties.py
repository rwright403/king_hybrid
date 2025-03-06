from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical
import matplotlib.pyplot as plt
import traceback


T_atm = 283 #K
P_atm = 1e5 #Pa

try:
    air = Chemical('air', T=T_atm, P=P_atm)
    k_f = air.kg
    visc_f = air.mug
    print("conductivity and viscosity of air from chemical lib: ", k_f, visc_f)
except Exception as e:
    print("thermo failed!")

import CoolProp.CoolProp as CP

T = 300  # Temperature in K
P = 101325  # Pressure in Pa

try:
    k_air = CP.PropsSI('L', 'T', T_atm, 'P', P_atm, 'Air')  # Thermal conductivity [W/m-K]
    mu_air = CP.PropsSI('V', 'T', T_atm, 'P', P_atm, 'Air')  # Dynamic viscosity [Pa.s]

    print(f"Thermal Conductivity of Air: {k_air:.6f} W/m-K")
    print(f"Viscosity of Air: {mu_air:.6e} Pa.s")

except Exception as e:
    print("coolprop failed!")