"""
import CoolProp.CoolProp as CP


T = 17 +273.15
rho_liq = CP.PropsSI('D', 'T', T, 'Q', 0.5, "N2O")
print(f"rho_liq = {rho_liq} kg/m^3")
"""

import CoolProp.CoolProp as CP
import numpy as np

# -----------------------
# Inputs
# -----------------------
P_tank = 6.5e6        # Tank pressure [Pa]
P_atm  = 101325       # Ambient pressure [Pa]

d_orifice = 0.0254*0.364    # Orifice diameter [m]
Cd = 1              # Discharge coefficient

fluid = 'NitrousOxide'

# -----------------------
# Geometry
# -----------------------
A = np.pi * (d_orifice/2)**2

# -----------------------
# Saturated vapor state (Q = 1)
# -----------------------
T_tank = CP.PropsSI('T', 'P', P_tank, 'Q', 1, fluid)

rho = CP.PropsSI('D', 'P', P_tank, 'Q', 1, fluid)   # density [kg/m^3]
a   = CP.PropsSI('A', 'P', P_tank, 'Q', 1, fluid)   # speed of sound [m/s]

# -----------------------
# Choked mass flow (sonic condition)
# -----------------------
mdot = Cd * A * rho * a

# -----------------------
# Exit conditions
# -----------------------
v_exit = a              # sonic velocity
P_exit = CP.PropsSI('P', 'T', T_tank, 'D', rho, fluid)  # ~tank pressure (approx)

# -----------------------
# Thrust
# -----------------------
F = mdot * v_exit + (P_exit - P_atm) * A

# -----------------------
# Output
# -----------------------
print(f"Saturation temperature: {T_tank:.2f} K")
print(f"Density: {rho:.2f} kg/m^3")
print(f"Speed of sound: {a:.2f} m/s")
print(f"Mass flow rate: {mdot:.3f} kg/s")
print(f"Thrust: {F:.2f} N")