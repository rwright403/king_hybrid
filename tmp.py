import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# --- Inputs ---
m_fuel = 50.0           # kg of ethanol
rho_fuel = 789.0        # kg/m^3
ullage_fraction = 0.15  # 15% ullage
fluid_gas = "Nitrogen"  # ullage gas

# --- Compute fuel and ullage volumes ---
V_fuel = m_fuel / rho_fuel
V_tank = V_fuel / (1 - ullage_fraction)
V_ullage = V_tank - V_fuel

# --- Gas properties from CoolProp ---
R_univ = 8.314462618  # J/(mol·K)
M_molar = PropsSI('M', fluid_gas)  # kg/mol
R_specific = R_univ / M_molar      # J/(kg·K)

# --- Temperature range ---
T_range = np.linspace(273, 323, 100)  # K, e.g., 0°C to 50°C

# --- Assume a fixed mass of nitrogen ---
# You can pick a reference pressure at a reference temperature to get mass
P_ref = 2e5  # Pa, reference pressure
T_ref = 293.15  # K, reference temperature
m_nitrogen = P_ref * V_ullage / (R_specific * T_ref)

# --- Compute pressure at each temperature ---
P = m_nitrogen * R_specific * T_range / V_ullage  # Pa

# --- Plot ---
plt.figure(figsize=(8,5))
plt.plot(T_range-273.15, P/1e5, color='blue')  # x in °C, y in bar
plt.xlabel("Temperature [°C]")
plt.ylabel("Nitrogen Pressure [bar]")
plt.title("Ullage Gas Pressure vs Temperature")
plt.grid(True)
plt.show()
