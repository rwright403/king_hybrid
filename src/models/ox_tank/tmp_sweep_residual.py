import numpy as np
import matplotlib.pyplot as plt
from src.models._thermo.n2o_thermo_span_wagner_n2o_constants import *

def lightweight_span_wagner_eos_pressure(rho, T): #NOTE: RETURNS P IN [Pa]
    tau = T_c / T
    delta = rho / rho_c
    ar_delta = np.sum(n1 * d1 * delta**(d1 - 1) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    return rho * R_mass * T * (1 + delta * ar_delta)

# Example tank state (replace with real values)
m_liq, m_gas = 18.0, 2.0     # kg
V_tank = 0.0354              # m^3
T_liq = 286.5                # K
T_gas = 286.5                # K
P_tank = 4.5e6               # Pa (clamped to Psat for demo)

def residual_rho_g(rho_g):
    V_gas = m_gas / rho_g
    V_liq = V_tank - V_gas
    if V_liq <= 0:
        return np.nan
    rho_liq = m_liq / V_liq
    P_liq = lightweight_span_wagner_eos_pressure(rho_liq, T_liq)
    return P_liq - P_tank

# Sweep over rho_g
rho_vals = np.linspace(0.1, 205, 500)
residuals = [residual_rho_g(r) for r in rho_vals]

# Plot
plt.figure(figsize=(8,5))
plt.axhline(0, color='k', linestyle='--', label="Residual=0")
plt.plot(rho_vals, residuals, label="P_liq - P_tank")
plt.xlabel(r"Gas density $\rho_g$ [kg/m³]")
plt.ylabel(r"Residual (Pa)")
plt.title("Residual sweep for P_liq - P_tank vs rho_g")
plt.legend()
plt.grid(True)
plt.show()


"""


def residual_rho_g(rho_g):
    V_gas = m_gas / rho_g
    V_liq = V_tank - V_gas
    rho_liq = m_liq / V_liq

    P_liq = lightweight_span_wagner_eos_pressure(rho_liq, T_liq)

    print("pressures: ", P_liq - P_tank, P_liq, P_tank)
    return P_liq - P_tank"""


