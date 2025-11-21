import numpy as np
from src.models.inj._base import BaseInjector
from src.models.inj.spi import spi_model
from src.models.inj.hem import hem_model

class dyer_model(BaseInjector):

    def __init__(self, Cd: float, A_inj: float):
        super().__init__(Cd, A_inj)
        self.spi = spi_model(Cd, A_inj)
        self.hem = hem_model(Cd, A_inj)


    def m_dot(self, state: dict) -> float:

        P_1 = state["P_1"]
        P_2 = state["P_2"]
        P_sat = state["P_sat"]

        # First compute SPI
        m_dot_spi = self.spi.m_dot(state)
                        
        #NOTE: FOR THIS CASE NO CAVITATION SO WE ARE JUST USING THE SPI MODEL
        if(P_sat < P_2):
            m_dot_dyer = m_dot_spi

        #NOTE: ELSE TWO PHASE AT INJ OUTLET AND USE DYER TO ACCOUNT FOR TWO PHASE EFFECTS
        else:
            m_dot_hem = self.hem.m_dot(state)
            
            dyer_k = np.sqrt( (P_1 - P_2) / (P_sat - P_2) ) 
            m_dot_dyer = ((dyer_k/(1+dyer_k)) * m_dot_spi) + ((1/(1+dyer_k)) * m_dot_hem)
        
        return m_dot_dyer

    


# === Quick test harness for Dyer model (pressure-based) ===
if __name__ == "__main__":
    import CoolProp.CoolProp as CP
    import numpy as np
    import matplotlib.pyplot as plt

    # Injector parameters
    Cd = 0.9
    d_inj = 1.5e-3   # [m]
    A_inj = np.pi * (d_inj / 2)**2

    # Initialize model
    inj = dyer_model(Cd, A_inj)

    # Upstream state (tank): define by pressure instead of temperature
    fluid = "N2O"
    P_1 = 3e6  # [Pa] — example tank pressure (50 bar)
    T_1 = CP.PropsSI("T", "P", P_1, "Q", 0, fluid)   # saturated liquid temperature
    rho_1 = CP.PropsSI("D", "P", P_1, "Q", 0, fluid)
    h_1 = CP.PropsSI("H", "P", P_1, "Q", 0, fluid)
    u_1 = CP.PropsSI("U", "P", P_1, "Q", 0, fluid)
    s_1 = CP.PropsSI("S", "P", P_1, "Q", 0, fluid)
    P_sat = P_1  # same since we’re on saturation line

    # Sweep downstream pressures
    P_2_arr = np.linspace(P_1 * 0.99, P_1 * 0.05, 80)
    m_dot_arr = np.zeros_like(P_2_arr)

    # Loop through downstream pressures
    for i, P_2 in enumerate(P_2_arr):
        state = {
            "P_1": P_1,
            "P_2": P_2,
            "T_1": T_1,
            "rho_1": rho_1,
            "h_1": h_1,
            "u_1": u_1,
            "s_1": s_1,
            "P_sat": P_sat,
            "x_1": 0.0,  # pure liquid inlet
        }
        try:
            m_dot_arr[i] = inj.m_dot(state)
        except Exception as e:
            print(f"Error at P₂ = {P_2/1e5:.2f} bar: {e}")
            m_dot_arr[i] = np.nan

    # Convert to convenient units
    deltaP_arr = (P_1 - P_2_arr) / 1e5  # [bar]
    m_dot_gs = m_dot_arr * 1e3           # [g/s]

    # === Plot: ṁ vs ΔP ===
    plt.figure(figsize=(6,4))
    plt.plot(deltaP_arr, m_dot_gs, "o-", lw=1.5)
    plt.xlabel("ΔP = P₁ - P₂  [bar]")
    plt.ylabel("ṁ  [g/s]")
    plt.title("Dyer Injector Model — Mass Flow vs ΔP")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # === Plot: ṁ vs P₂ to visualize choking ===
    plt.figure(figsize=(6,4))
    plt.plot(P_2_arr / 1e5, m_dot_gs, "o-", lw=1.5)
    plt.xlabel("Downstream Pressure P₂ [bar]")
    plt.ylabel("ṁ [g/s]")
    plt.title("Dyer Model — Choked Flow Behavior")
    plt.gca().invert_xaxis()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
