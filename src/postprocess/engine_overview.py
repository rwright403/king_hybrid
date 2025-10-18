import numpy as np

def engine_overview(program_input, results, mode):

    if mode == "full_stack":

        g0 = 9.80665  # gravitational acceleration [m/s^2]

        time = results["time"]
        thrust = results["thrust"]
        m_dot_cc = results["m_dot_cc"]

        peak_thrust = np.max(thrust)

        # eqn (2-3) from rpe: take integral of thrust wrt time / (g_0 * integral of nozzle masss flow wrt time)

        # Compute total impulse: integral of thrust over time
        total_impulse = np.trapezoid(thrust, time)  # [N·s]
        

        # Compute total propellant mass expelled
        total_mass_expended = np.trapezoid(m_dot_cc, time)  # [kg]

        # Compute specific impulse
        Isp = total_impulse / (g0 * total_mass_expended)

        print(f"\n\n--- ENGINE PERFORMANCE OVERVIEW ---")
        print(f"Total Impulse (I_t):        {total_impulse:.2f} N·s")
        print(f"Total Mass Expended:        {total_mass_expended:.4f} kg")
        print(f"Specific Impulse (I_sp):    {Isp:.2f} s")
        print(f"Peak Thrust:                {peak_thrust:.2f} N")
        print(f"Burn Time:                  {program_input.sim_time:.2f} s, (simulated time)")


# test harness:



# -----------------------
# Test harness section
# -----------------------
"""
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class ProgramInput:
    sim_time: float   # [s]


def run_engine_overview_from_split_csvs(thrust_csv, mdot_csv):
    #Reads separate thrust and m_dot_cc CSVs, merges into a single results dict, and runs engine_overview().

    # Read both CSVs, skipping header row (row 0)
    thrust_df = pd.read_csv(thrust_csv, skiprows=1, names=["time", "thrust"])
    mdot_df   = pd.read_csv(mdot_csv,   skiprows=1, names=["time", "m_dot_cc"])

    # Validate equal lengths and matching time
    if len(thrust_df) != len(mdot_df):
        raise ValueError("Thrust and mass flow CSVs have different lengths!")
    if not np.allclose(thrust_df["time"], mdot_df["time"], atol=1e-9):
        raise ValueError("Time columns in thrust and mass flow CSVs do not match!")

    # Merge into a single DataFrame
    df = pd.DataFrame({
        "time": thrust_df["time"],
        "thrust": thrust_df["thrust"],
        "m_dot_cc": mdot_df["m_dot_cc"],
    })

    # Optional: plot for debugging

    fig, ax1 = plt.subplots()
    ax1.plot(df["time"], df["thrust"], label="Thrust [N]", linewidth=2)
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Thrust [N]", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")

    # second y-axis for mass-flow
    ax2 = ax1.twinx()
    ax2.plot(df["time"], df["m_dot_cc"], color="tab:orange", linestyle="--", label="Mass Flow [kg/s]")
    ax2.set_ylabel("Mass Flow [kg/s]", color="tab:orange")
    ax2.tick_params(axis="y", labelcolor="tab:orange")

    plt.title("Thrust & Mass Flow vs Time")
    fig.tight_layout()
    plt.show()

    return df



# Replace with actual file paths:
thrust_csv_path = "src/results/rocketpy_test_UofT_style_rocket/thrust.csv"
mdot_csv_path = "src/results/rocketpy_test_UofT_style_rocket/m_dot_cc.csv"

df_combined = run_engine_overview_from_split_csvs(thrust_csv_path, mdot_csv_path)

program_input = ProgramInput(sim_time=5.0)

engine_overview(program_input, results=df_combined, mode="full_stack")
"""