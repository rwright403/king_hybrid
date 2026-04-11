import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def load_and_interp(csv_path, t_common, value_col="value"):
    df = pd.read_csv(csv_path)
    t = df["time"].values
    y = df[value_col].values
    interp = interp1d(
        t, y,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate"
    )
    return interp(t_common)

def nozzle_performance(
    ve_csv,
    mdot_csv,
    pc_csv,
    pe_csv,
    pa_csv,
    alt_csv,
    At,
    Ae,
    t_burn,
    n_points=2000
):
    """
    Compute and plot Cf vs time over burn duration only.
    """
    # --- Load CSVs ---
    t_start = 0.0
    t_end = t_burn
    t_common = np.linspace(t_start, t_end, n_points)

    ve = load_and_interp(ve_csv, t_common, value_col="value")
    mdot = load_and_interp(mdot_csv, t_common, value_col="value")
    pc = load_and_interp(pc_csv, t_common, value_col="value")
    pe = load_and_interp(pe_csv, t_common, value_col="value")
    pa = load_and_interp(pa_csv, t_common, value_col="value")
    alt = load_and_interp(alt_csv, t_common, value_col="value")

    # --- Compute thrust coefficient ---
    delta_p = (pe - pa)
    thrust = mdot * ve + (delta_p) * Ae
    Cf = thrust / (pc * At)


    # --- Plot ---
    plt.subplot(1, 3, 1)
    plt.plot(t_common, Cf)
    plt.xlabel("Time [s]")
    plt.ylabel("Thrust Coefficient $C_f$")
    plt.title("Nozzle Thrust Coefficient During Burn")
    plt.grid(True)

    plt.subplot(1, 3, 2)
    plt.plot(alt, Cf)
    plt.xlabel("Alt [m]")
    plt.ylabel("Thrust Coefficient $C_f$")
    plt.title("Nozzle Thrust Coefficient vs Altitude")
    plt.grid(True)

    plt.subplot(1, 3, 3)
    plt.plot(alt, delta_p)
    plt.xlabel("Alt [m]")
    plt.ylabel("Delta_P [Pa]")
    plt.title("P_exit - P_atm vs Altitude")
    plt.grid(True)

    plt.show()