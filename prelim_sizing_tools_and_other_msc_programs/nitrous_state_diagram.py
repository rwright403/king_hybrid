import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

# Define the substance
substance = 'NitrousOxide'

# Create a range of temperatures for the saturated region
T_min, T_max = CP.PropsSI('Tmin', substance), CP.PropsSI('Tcrit', substance)
temperatures = np.linspace(T_min, T_max, 500)

# Initialize arrays for saturated properties
P_sat = np.zeros_like(temperatures)
V_sat_liq = np.zeros_like(temperatures)
V_sat_vap = np.zeros_like(temperatures)
S_sat_liq = np.zeros_like(temperatures)
S_sat_vap = np.zeros_like(temperatures)

# Calculate saturated properties
for i, T in enumerate(temperatures):
    P_sat[i] = CP.PropsSI('P', 'T', T, 'Q', 0, substance)
    V_sat_liq[i] = 1 / CP.PropsSI('D', 'T', T, 'Q', 0, substance)
    V_sat_vap[i] = 1 / CP.PropsSI('D', 'T', T, 'Q', 1, substance)
    S_sat_liq[i] = CP.PropsSI('S', 'T', T, 'Q', 0, substance) / 1000  # Convert J/(kg*K) to kJ/(kg*K)
    S_sat_vap[i] = CP.PropsSI('S', 'T', T, 'Q', 1, substance) / 1000  # Convert J/(kg*K) to kJ/(kg*K)

# Create subplots for P-V and T-S diagrams side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))

# Plot P-V diagram
ax1.plot(V_sat_liq, P_sat, 'b-', label='Saturated Liquid')
ax1.plot(V_sat_vap, P_sat, 'r-', label='Saturated Vapor')
ax1.set_title('Pressure-Volume Diagram for N2O')
ax1.set_xlabel('Volume (m^3/kg)')
ax1.set_ylabel('Pressure (Pa)')
#ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend(loc='best')
ax1.grid(True, which="both", ls="--")

#ax1.axhline(y=CP.PropsSI('P', 'T', 282, 'Q', 1, substance), color='rebeccapurple', linestyle='--')
ax1.axhline(y=(2.55e6 + CP.PropsSI('P', 'T', 282, 'Q', 1, substance)), color='rebeccapurple', linestyle='--')
#288.1145796254294 4480926.129652254
ax1.plot( 1/CP.PropsSI('D', 'T', 288.1145796254294, 'P', 45e5+100, substance), 45e5+100, marker="o", markersize=10, color="red")
ax1.plot( 1/CP.PropsSI('D', 'T', 288.1145796254294, 'P', 4480926.129652254, substance), 4480926.129652254, marker="o", markersize=10, color="rebeccapurple")
#ax1.scatter(y = (2.55e6 + CP.PropsSI('P', 'T', T, 'Q', 282, 'N2O') ), )

# Plot T-S diagram
ax2.plot(S_sat_liq, temperatures, 'b-', label='Saturated Liquid')
ax2.plot(S_sat_vap, temperatures, 'r-', label='Saturated Vapor')
ax2.set_title('Temperature-Entropy Diagram for N2O')
ax2.set_xlabel('Entropy (kJ/(kg*K))')
ax2.set_ylabel('Temperature (K)')
ax2.legend(loc='best')
ax2.grid(True, which="both", ls="--")

plt.tight_layout()
plt.show()
