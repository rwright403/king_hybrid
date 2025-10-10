import numpy as np
import CoolProp.CoolProp as CP
from src.models.inj._base import BaseInjector
from src.models.inj.spi import spi_model
from src.models.inj.hem import hem_model


class dyer_model(BaseInjector):

    def __init__(self, Cd: float, A_inj: float):
        super().__init__(Cd, A_inj)
        self.spi = spi_model(Cd, A_inj)
        self.hem = hem_model(Cd, A_inj)

    def m_dot(self, state: dict) -> float:
        P_1  = state["P_1"]
        P_2  = state["P_2"]
        P_sat = state["P_sat"]
        h_1  = state["h_1"]

        # --- SPI baseline
        m_dot_spi = self.spi.m_dot(state)

        # --- Single-phase regime
        if P_sat < P_2:
            return m_dot_spi

        # --- Two-phase regime → compute Dyer blend
        m_dot_hem = self.hem.m_dot(state)
        dyer_k = np.sqrt(max(P_1 - P_2, 0.0) / max(P_sat - P_2, 1e-6))
        m_dot_dyer = (dyer_k / (1.0 + dyer_k)) * m_dot_spi + (1.0 / (1.0 + dyer_k)) * m_dot_hem

        # ============================================================
        # HEM-BASED CHOKING CHECK (internal)
        # ============================================================

        # upstream entropy (isentropic assumption)
        s_1 = CP.PropsSI("S", "H", h_1, "P", P_1, "N2O")

        # build isentropic expansion curve from P1→P2
        downstream_pres_arr = np.linspace(P_2, P_1, 100)
        h_2_list = []
        rho_2_list = []

        for pres in downstream_pres_arr:
            h_2_list.append(CP.PropsSI("H", "S", s_1, "P", pres, "N2O"))
            rho_2_list.append(CP.PropsSI("D", "S", s_1, "P", pres, "N2O"))

        h_2_arr = np.array(h_2_list)
        rho_2_arr = np.array(rho_2_list)

        # compute corresponding HEM mass flow curve
        m_dot_arr = self.C_inj * rho_2_arr * np.sqrt(2.0 * np.abs(h_1 - h_2_arr))

        # find the critical (choked) condition
        idx_crit = np.argmax(m_dot_arr)
        m_dot_crit = m_dot_arr[idx_crit]
        P_crit = downstream_pres_arr[idx_crit]

        # if downstream pressure is below the critical point → choked
        if P_2 <= P_crit:
            m_dot_dyer = min(m_dot_dyer, m_dot_crit)

        return m_dot_dyer

    


# === Quick test harness for Dyer model ===
if __name__ == "__main__":
    import CoolProp.CoolProp as CP
    import numpy as np
    import matplotlib.pyplot as plt

    # Injector parameters
    Cd = 0.9
    d_inj = 1.5e-3   # [m]
    A_inj = np.pi * (d_inj/2)**2

    # Initialize model
    inj = dyer_model(Cd, A_inj)

    # Upstream (tank) state: choose a nominal liquid N2O condition
    T_1 = 293.15  # [K]
    fluid = "N2O"
    P_sat = CP.PropsSI("P", "T", T_1, "Q", 0, fluid)
    rho_1 = CP.PropsSI("D", "T", T_1, "Q", 0, fluid)
    h_1 = CP.PropsSI("H", "T", T_1, "Q", 0, fluid)
    s_1 = CP.PropsSI("S", "T", T_1, "Q", 0, fluid)
    u_1 = CP.PropsSI("U", "T", T_1, "Q", 0, fluid)

    # Slightly pressurized above saturation (subcooled liquid)
    P_1 = 1.05 * P_sat

    # Sweep downstream pressures
    P_2_arr = np.linspace(P_1*0.99, P_1*0.05, 80)
    m_dot_arr = np.zeros_like(P_2_arr)

    # Loop through pressures
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
            "x_1": 0.0,
        }
        try:
            m_dot_arr[i] = inj.m_dot(state)
        except Exception as e:
            print(f"Error at P2={P_2/1e5:.2f} bar: {e}")
            m_dot_arr[i] = np.nan

    # Convert to convenient units
    deltaP_arr = (P_1 - P_2_arr) / 1e5  # [bar]
    m_dot_gs = m_dot_arr * 1e3          # [g/s]

    # === Plot results ===
    plt.figure(figsize=(6,4))
    plt.plot(deltaP_arr, m_dot_gs, "o-", lw=1.5)
    plt.xlabel("ΔP = P₁ - P₂  [bar]")
    plt.ylabel("ṁ  [g/s]")
    plt.title("Dyer Injector Model — Mass Flow vs ΔP")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Optionally plot ṁ vs P₂ to visualize choke plateau
    plt.figure(figsize=(6,4))
    plt.plot(P_2_arr/1e5, m_dot_gs, "o-", lw=1.5)
    plt.xlabel("Downstream Pressure P₂ [bar]")
    plt.ylabel("ṁ [g/s]")
    plt.title("Dyer Model — Choked Flow Behavior")
    plt.gca().invert_xaxis()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
