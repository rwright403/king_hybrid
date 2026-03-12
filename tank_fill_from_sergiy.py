"""
N2O Fill Time Surface Plot – Two-Phase Tank with Energy Balance (Parallel, Choked Check)

Features:
    • Rigid N2O tank with mass + energy balance (flash, chilling).
    • Real-fluid thermo via CoolProp (NitrousOxide).
    • Uniform two-phase control volume (liquid + vapour).
    • Vent modeled as compressible gas orifice, with choked vs unchoked check.
    • 3D surface: fill time vs vent orifice diameter vs inlet mass flow.
    • Parallel grid computation using ProcessPoolExecutor.
    • tqdm progress bar for total progress.

Limitations:
    • Vent flow treated as vapour (top ullage vent assumption), no explicit two-phase choking.
    • No line losses, no wall heat transfer, no injector-side modeling.

Dependencies (inside your Python 3.11 venv):
    pip install numpy matplotlib CoolProp tqdm
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from matplotlib.colors import LogNorm
from CoolProp.CoolProp import PropsSI
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# ---------------------------------------------------------
# 1. GLOBAL PARAMETERS & CONSTANTS
# ---------------------------------------------------------

FLUID = "NitrousOxide"  # CoolProp fluid name for N2O

# Tank geometry
V_TANK = 0.031467        # [m^3] internal tank volume (given)
FILL_FRACTION_TARGET = 0.75  # [-] target liquid volume fraction (75%)

# Liquid N2O density at ~20 °C [kg/m^3]
# Airgas SDS gives ~0.785 kg/L at 20 °C => 785 kg/m^3:
# https://www.airgasspecialtyproducts.com/wp-content/uploads/2016/04/Nitrous-Oxide-SDS.pdf
RHO_LIQ_REF = 785.0

DELTA_M_TARGET = RHO_LIQ_REF * (FILL_FRACTION_TARGET * V_TANK)  # [kg]

# N2O vapour properties (approx constants near 20 °C)
# cp ≈ 0.88 kJ/kg/K, cv ≈ 0.69 kJ/kg/K:
# https://www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html
CP_KJ = 0.88
CV_KJ = 0.69
GAMMA = CP_KJ / CV_KJ                 # [-] cp/cv
R_SPEC = (CP_KJ - CV_KJ) * 1000.0    # [J/(kg·K)]

# Reference saturation state for initial condition & inlet
T_REF = 293.15                        # [K] ~20 °C
# Vapour pressure at 20 °C ≈ 50.6 bar:
# https://www.engineeringtoolbox.com/nitrous-oxide-d_2034.html
P_REF = 50.6e5                        # [Pa]

# Ambient / receiver pressure (vent to atmosphere)
P_AMB = 1.0e5                         # [Pa]

# Orifice discharge coefficient
C_D = 0.8                             # [-]

# Time-stepping parameters
DT = 0.05                             # [s]
T_MAX = 600.0                         # [s] max simulated time


# ---------------------------------------------------------
# 2. PARAMETER RANGES (VENT DIAMETER & INLET MASS FLOW)
# ---------------------------------------------------------

D_MIN = 0.002                         # [m] 2 mm
D_MAX = 0.010                         # [m] 10 mm
N_D = 20							  #RESOLUTION CONST DIAMETER

M_DOT_IN_MIN = 0.1                    # [kg/s]
M_DOT_IN_MAX = 5.0                    # [kg/s]
N_M_DOT = 20						  #RESOLUTION CONST MASS FLOW


# ---------------------------------------------------------
# 3. HELPER FUNCTIONS – THERMODYNAMICS
# ---------------------------------------------------------

def initial_tank_state():
    """
    Initial state: saturated vapour at T_REF.
    """
    rho_v = PropsSI("Dmass", "T", T_REF, "Q", 1.0, FLUID)   # [kg/m^3]
    u_v = PropsSI("Umass", "T", T_REF, "Q", 1.0, FLUID)     # [J/kg]
    m0 = rho_v * V_TANK                                     # [kg]
    U0 = m0 * u_v                                           # [J]
    return m0, U0


def state_from_mU(m, U):
    """
    Given total mass m [kg] and internal energy U [J], get P, T, Q.

    Uses:
        rho = m / V_TANK
        u   = U / m
    and CoolProp with (Dmass, Umass) inputs.
    """
    rho = m / V_TANK                          # [kg/m^3]
    u = U / m                                 # [J/kg]

    P = PropsSI("P", "Dmass", rho, "Umass", u, FLUID)  # [Pa]
    T = PropsSI("T", "Dmass", rho, "Umass", u, FLUID)  # [K]
    Q = PropsSI("Q", "Dmass", rho, "Umass", u, FLUID)  # [-]

    return P, T, Q


def liquid_volume_fraction(m, P, T, Q):
    """
    Compute liquid volume fraction in tank.

    If 0 < Q < 1: two-phase mixture. Otherwise treat as fully liquid / fully vapour.
    """
    if Q <= 0.0:
        return 1.0
    if Q >= 1.0:
        return 0.0

    rho_l = PropsSI("Dmass", "T", T, "Q", 0.0, FLUID)  # [kg/m^3]
    rho_v = PropsSI("Dmass", "T", T, "Q", 1.0, FLUID)  # [kg/m^3]

    m_v = Q * m
    m_l = (1.0 - Q) * m

    V_l = m_l / rho_l
    V_v = m_v / rho_v

    return V_l / (V_l + V_v)


# ---------------------------------------------------------
# 4. VENT MASS FLOW (CHOKED VS UNCHOKED CHECK)
# ---------------------------------------------------------

def vent_mass_flow(diameter_m, P_up, T_up):
    """
    Mass flow through vent orifice, with choked/un-choked decision.

    Compressible (ideal-gas) orifice relations:

        p_crit  = (2/(γ+1))^(γ/(γ−1))
        p_ratio = P_down / P_up

    If p_ratio <= p_crit (choked):
        m_dot = C_d * A * P_up / sqrt(T_up*R) *
                sqrt(γ) * (2/(γ+1))^((γ+1)/(2(γ−1)))

    Else (subcritical):
        m_dot = C_d * A * P_up / sqrt(T_up*R) *
                sqrt( 2γ/(γ−1) * [(p_ratio)^(2/γ) − (p_ratio)^((γ+1)/γ)] )

    P_down = P_AMB = 1 bar (vent to atmosphere).
    """

    A = 0.25 * np.pi * diameter_m**2      # [m^2] orifice area
    P_down = P_AMB                        # [Pa]

    p_ratio = P_down / P_up
    p_crit = (2.0 / (GAMMA + 1.0)) ** (GAMMA / (GAMMA - 1.0))

    prefactor = C_D * A * P_up / np.sqrt(T_up * R_SPEC)

    if p_ratio <= p_crit:
        # Choked flow
        factor = np.sqrt(GAMMA) * (2.0 / (GAMMA + 1.0)) ** ((GAMMA + 1.0) /
                                                            (2.0 * (GAMMA - 1.0)))
        m_dot = prefactor * factor
    else:
        # Subcritical (unchoked)
        bracket = (2.0 * GAMMA / (GAMMA - 1.0) *
                   (p_ratio ** (2.0 / GAMMA) -
                    p_ratio ** ((GAMMA + 1.0) / GAMMA)))
        bracket = max(bracket, 0.0)
        m_dot = prefactor * np.sqrt(bracket)

    return m_dot


# ---------------------------------------------------------
# 5. FILL SIMULATION FOR ONE (DIAMETER, M_DOT_IN) PAIR
# ---------------------------------------------------------

# Inlet enthalpy: saturated liquid at T_REF
H_IN = PropsSI("Hmass", "T", T_REF, "Q", 0.0, FLUID)   # [J/kg]


def simulate_fill_time(diameter_m, m_dot_in,
                       dt=DT, t_max=T_MAX):
    """
    Simulate one fill transient for a given vent diameter and inlet mass flow.

    Mass balance:
        dm/dt = m_dot_in - m_dot_out

    Energy balance:
        dU/dt = m_dot_in * h_in - m_dot_out * h_out

    Stop when liquid volume fraction >= FILL_FRACTION_TARGET.
    Return NaN if net flow is negative long-term or t_max exceeded.
    """
    m, U = initial_tank_state()
    t = 0.0

    m_prev = m

    while t < t_max:
        if m <= 1e-8:
            return np.nan

        try:
            P, T, Q = state_from_mU(m, U)
        except Exception:
            return np.nan

        phi_l = liquid_volume_fraction(m, P, T, Q)

        if phi_l >= FILL_FRACTION_TARGET:
            return t

        # Vent mass flow
        m_dot_out = vent_mass_flow(diameter_m, P, T)

        m_dot_net = m_dot_in - m_dot_out

        # If we're clearly blowing down rather than filling
        if (t > 5.0) and (m < m_prev) and (m_dot_net <= 0.0):
            return np.nan

        # Enthalpy of vented saturated vapour at current P
        h_out = PropsSI("Hmass", "P", P, "Q", 1.0, FLUID)  # [J/kg]

        dm = m_dot_net * dt
        dU = (m_dot_in * H_IN - m_dot_out * h_out) * dt

        m_prev = m
        m += dm
        U += dU
        t += dt

        if m <= 0.0:
            return np.nan

    return np.nan


# ---------------------------------------------------------
# 6. PARALLEL GRID EVALUATION + PLOTTING (MAIN)
# ---------------------------------------------------------

def _worker(job):
    """
    Worker function for ProcessPoolExecutor.

    job = (i, j, d, m_in)
    """
    i, j, d, m_in = job
    t_fill = simulate_fill_time(d, m_in)
    return i, j, t_fill


def main():
    # Build parameter arrays
    d_vals = np.linspace(D_MIN, D_MAX, N_D)                        # [m]
    m_in_vals = np.linspace(M_DOT_IN_MIN, M_DOT_IN_MAX, N_M_DOT)   # [kg/s]

    D_mesh, M_in_mesh = np.meshgrid(d_vals, m_in_vals)
    T_fill = np.empty_like(D_mesh, dtype=float)

    # Build job list of (i, j, d, m_in)
    jobs = []
    for i in range(N_M_DOT):
        for j in range(N_D):
            jobs.append((i, j, d_vals[j], m_in_vals[i]))

    print("---------------------------------------------------")
    print("N₂O Two-Phase Fill Model Starting (Parallel)…")
    print(f"Grid size: {N_M_DOT} inlet-flow points × {N_D} vent diameters")
    print(f"Total simulations: {len(jobs)}")
    print(f"Using up to {multiprocessing.cpu_count()} CPU cores.\n")
    print("---------------------------------------------------")

    # Parallel execution with progress bar
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as ex:
        futures = [ex.submit(_worker, job) for job in jobs]
        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc="Simulating grid", unit="case"):
            i, j, t_val = fut.result()
            T_fill[i, j] = t_val

    print("\n---------------------------------------------------")
    print("Surface computation complete.")
    print("---------------------------------------------------\n")
    
    # Mask invalid / non-positive fill times
    T_plot = np.where(T_fill > 0.0, T_fill, np.nan)

    # -------------------------------------------------
    # 7a. 3D SURFACE PLOT (same as before)
    # -------------------------------------------------
    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection="3d")

    ax.set_zscale("log")

    surf = ax.plot_surface(
        D_mesh * 1000.0,      # [mm]
        M_in_mesh,            # [kg/s]
        T_plot,               # [s]
        cmap="viridis",
        linewidth=0,
        antialiased=True
    )

    z_min_finite = np.nanmin(T_plot)
    ax.contour(
        D_mesh * 1000.0,
        M_in_mesh,
        T_plot,
        levels=15,
        cmap="viridis",
        offset=z_min_finite
    )

    ax.set_xlabel("Vent orifice diameter [mm]")
    ax.set_ylabel("Fill mass flow rate ṁ_in [kg/s]")
    ax.set_zlabel("Fill time to 75% liquid volume [s] (log scale)")
    ax.set_title("N₂O Tank Fill Time vs Vent Orifice Diameter and Fill Mass Flow\n"
                 "(Two-phase tank + energy balance, choked/unchoked vent model)")

    cbar = fig.colorbar(surf, shrink=0.55, aspect=14)
    cbar.set_label("Fill time [s] (log scale)")

    # -------------------------------------------------
    # 7b. 2D CONTOUR MAP (top-down design map)
    # -------------------------------------------------
    fig2, ax2 = plt.subplots(figsize=(8, 6))

    # Use log-spaced contour levels based on valid times
    valid = np.isfinite(T_plot)
    if np.any(valid):
        t_min = np.nanmin(T_plot[valid])
        t_max = np.nanmax(T_plot[valid])
        # Guard against degenerate ranges
        if t_max <= t_min * 1.01:
            t_max = t_min * 1.01

        levels = np.logspace(np.floor(np.log10(t_min)),
                             np.ceil(np.log10(t_max)),
                             15)

        cf = ax2.contourf(
            D_mesh * 1000.0,
            M_in_mesh,
            T_plot,
            levels=levels,
            norm=LogNorm(),
            cmap="viridis"
        )
        cs = ax2.contour(
            D_mesh * 1000.0,
            M_in_mesh,
            T_plot,
            levels=levels,
            norm=LogNorm(),
            colors="k",
            linewidths=0.5
        )

        cbar2 = fig2.colorbar(cf)
        cbar2.set_label("Fill time [s] (log scale)")

        ax2.clabel(cs, inline=True, fmt="%.0f s", fontsize=8)
    else:
        ax2.text(0.5, 0.5, "No valid data", ha="center", va="center",
                 transform=ax2.transAxes)

    ax2.set_xlabel("Vent orifice diameter [mm]")
    ax2.set_ylabel("Fill mass flow rate ṁ_in [kg/s]")
    ax2.set_title("N₂O Tank Fill Time – Contour Map\n"
                  "(log-scaled time levels)")

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------
# ENTRY POINT (REQUIRED FOR MULTIPROCESSING ON macOS)
# ---------------------------------------------------------

if __name__ == "__main__":
    main()
