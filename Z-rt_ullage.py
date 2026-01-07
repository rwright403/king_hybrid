import numpy as np
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt

def fill_fraction(T, V_tank, m_total, fluid="N2O"):
    rho_l = CP.PropsSI("D","T",T,"Q",0,fluid)
    rho_v = CP.PropsSI("D","T",T,"Q",1,fluid)

    # quality x
    x = (V_tank/m_total - 1/rho_l) / (1/rho_v - 1/rho_l)
    x = np.clip(x, 0, 1)

    # liquid volume
    V_liq = (1 - x) * m_total / rho_l

    # fill fraction
    return V_liq / V_tank


"""
def gas_height_equilibrium(T, V_tank, m_total, diameter):

    rho_l = CP.PropsSI("D", "T", T, "Q", 0, "N2O")  # saturated liquid density
    rho_v = CP.PropsSI("D", "T", T, "Q", 1, "N2O")  # saturated vapor density

    # Solve vapor mass quality x
    x = (V_tank/m_total - 1/rho_l) / (1/rho_v - 1/rho_l)

    # Masses
    m_v = x * m_total
    m_l = m_total - m_v

    # Volumes
    V_v = m_v / rho_v
    V_l = m_l / rho_l

    # Heights
    A = np.pi * (diameter/2)**2
    h_gas = V_v / A
    h_liq = V_l / A

    print(f"h_gas: {h_gas:.3f} m, h_liq: {h_liq:.3f} m")

    #return h_gas, h_liq
"""


m_ox = 7.650873122104811 #kg
T_min = 273.15 + 17
T_max = 273.15 + 33
V_selected = 0.013
diam = 0.0254*5.0

# Parameters
fluid = "N2O"
V_min = m_ox / CP.PropsSI("D","T",T_max,"Q",0,fluid)  # min tank volume
V_max = V_min * 1.5 #allow tank volumes up to 1.75 times the min volume for 100% sat liquid.
V_points = np.linspace(V_min, V_max, 100)
T_points = np.linspace(T_min, T_max, 100)

# Compute 2D FF array
FF = np.zeros((len(T_points), len(V_points)))
P = np.zeros_like(FF)

for i, T in enumerate(T_points):
    rho_l = CP.PropsSI("D","T",T,"Q",0,fluid)
    rho_v = CP.PropsSI("D","T",T,"Q",1,fluid)

    # solve for quality x for every tank volume
    x = (V_points/m_ox - 1/rho_l) / (1/rho_v - 1/rho_l)

    # compute true liquid & vapor volumes
    V_liq = (1-x)*m_ox/rho_l
    # V_vap = x*m_ox/rho_v   # not needed for fill fraction

    # true fill fraction
    FF[i,:] = V_liq / V_points

    # pressure at saturated T
    P[i,:] = CP.PropsSI("P","T",T,"Q",0,fluid)

# Plot as color map
fig, ax = plt.subplots(figsize=(8,5))
c = ax.pcolormesh(V_points, P/1e6, FF, cmap='plasma', shading='auto')  # Pressure in MPa
# Move colorbar to the very right of the figure
cbar = fig.colorbar(c, ax=ax, orientation='vertical', fraction=0.015, pad=0.02)
cbar.set_label("Fill Fraction (FF)")

ax.set_xlabel("Tank Volume [m^3]")
ax.set_ylabel("Pressure [MPa]")
ax.set_title("Nitrous Tank Volume vs Pressure Color Map (Fill Fraction)")

# Compute pressures for T_min and T_max
P_min = CP.PropsSI("P","T",T_min,"Q",0,fluid) / 1e6  # MPa
P_max = CP.PropsSI("P","T",T_max,"Q",0,fluid) / 1e6  # MPa

# Draw horizontal dashed lines
ax.hlines([P_min, P_max], V_min, V_max, colors='red', linestyles='dashed')

# Add labels at the right end of the lines
ax.text(V_max*1.01, P_min, f'T_min = {T_min-273.15:.1f} °C', va='center', color='red')
ax.text(V_max*1.01, P_max, f'T_max = {T_max-273.15:.1f} °C', va='center', color='red')

ax.vlines(V_selected, P.min()/1e6, P.max()/1e6,
          colors='red', linestyles='dashed', linewidth=2)

# Add a label next to the line
ax.text(V_selected, (P.max()/1e6)*1.01,
        f"V = {V_selected:.4f} m³",
        color='red', ha='center')

#print("\nUse the graph to estimate the tank volume. Note fill fraction = V_liq/V_tank.\nWe want a tank that isn't too big (long fill time)\nbut we don't want our fill fraction to be too high (overpressurization risk if temp increases)")

#gas_height_equilibrium(T_min, V_selected, m_ox, diam)

# Compute fill fraction for V_selected across all temperatures
FF_min = fill_fraction(T_min, V_selected, m_ox)
FF_max = fill_fraction(T_max, V_selected, m_ox)

print("FF at Tmin:", FF_min)
print("FF at Tmax:", FF_max)




plt.show()

