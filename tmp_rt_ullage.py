from CoolProp.CoolProp import PropsSI
import numpy as np

def gas_height_equilibrium(T, V_tank, m_total, diameter):
    """
    Solve vapor height in a vertical cylindrical tank of nitrous oxide at equilibrium.
    
    Inputs:
        P        : pressure [Pa]
        T        : temperature [K]
        V_tank   : total tank volume [m^3]
        m_total  : total mass of nitrous in tank [kg]
        diameter : inner tank diameter [m]
        
    Returns:
        h_gas : height of vapor column [m]
        h_liq : height of liquid column [m]
    """

    # --- Step 1: Saturated liquid and vapor densities at T ---
    rho_liq = PropsSI("D", "T", T, "Q", 0, "N2O")  # saturated liquid density
    rho_vap = PropsSI("D", "T", T, "Q", 1, "N2O")  # saturated vapor density

    # --- Step 2: Solve for vapor volume using equilibrium constraint ---
    # mass = rho_liq*V_liq + rho_vap*V_vap
    # with V_liq = V_tank - V_vap
    # => m = rho_liq*(V_tank - V_vap) + rho_vap*V_vap
    # Solve explicitly:
    V_vap = (m_total - rho_liq * V_tank) / (rho_vap - rho_liq)

    # Clip to physical limits
    V_vap = max(0.0, min(V_vap, V_tank))

    # --- Step 3: Convert vapor volume to height ---
    radius = diameter / 2
    A_cross = np.pi * radius**2

    h_gas = V_vap / A_cross
    h_liq = (V_tank - V_vap) / A_cross

    return h_gas, h_liq


# --------------------- Example Usage ---------------------
if __name__ == "__main__":       
    T = 273.15 + 5           # K
    V_tank = 0.031467   # m^3
    m_total = 17.7    # kg
    diameter = 0.0254*6 #m

    h_gas, h_liq = gas_height_equilibrium(T, V_tank, m_total, diameter)

    print(f"Gas height:    {h_gas:.3f} m")
    print(f"Liquid height: {h_liq:.3f} m")
