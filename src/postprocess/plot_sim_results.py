import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from itertools import cycle

# base colors per key (extend as needed)
BASE_COLORS = {
    "P_ox_tank": "tab:blue",
    "m_ox_tank": "tab:purple",
    "m_dot_ox_tank": "tab:pink",
    "P_fuel_tank": "tab:orange",
    "P_cc": "tab:green",
    "thrust": "tab:red",
    "T_ox_tank": "tab:brown",
}

# fallback color cycle
COLOR_CYCLE = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])


def lighter_color(color, factor=0.5):
    """
    Lighten the given color by blending with white.
    factor=0 keeps original, factor=1 gives white.
    """
    if color is None:
        return None
    c = mcolors.to_rgb(color)
    return tuple(1 - (1 - x) * (1 - factor) for x in c)


def load_validation_data(filepaths: dict) -> dict:
    """
    Load validation data CSVs into pandas DataFrames.
    Skips header row, assumes 2 cols: time + value.
    Uses dict key to label the signal.
    """
    data = {}
    for name, path in filepaths.items():
        df = pd.read_csv(path, header=None, skiprows=1, names=["time", "value"])
        data[name] = df
    return data


def plot_series(time, results_dict, validation_dict=None, title=None, ylabel=None):
    for key, values in results_dict.items():
        # pick base color, fallback if not defined
        base_color = BASE_COLORS.get(key, None)
        if base_color is None:
            base_color = next(COLOR_CYCLE)

        # plot sim
        plt.plot(time, values, label=f"sim: {key}", color=base_color)

        # plot validation
        if validation_dict and key in validation_dict:
            df = validation_dict[key]
            val_color = lighter_color(base_color)
            plt.plot(df["time"], df["value"],label=f"val: {key}", color=val_color)

    plt.title(title if title else "Results")
    plt.xlabel("Time [s]")
    if ylabel:
        plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True)


def plot_sim_results(inputs, results: dict, mode: str, save_path: str = None):
    time = results["time"]

    # Load validation if available
    validation_data = None
    if getattr(inputs, "validation_files", None):
        try:
            validation_data = load_validation_data(inputs.validation_files)
        except Exception as e:
            print(f"Warning: validation load failed -> {e}")

    plt.figure(figsize=(12, 6))

    if mode == "ox_tank":
        if inputs.ox_tank_model == 1:  # equilibrium
            plt.subplot(1, 3, 1)
            plot_series(time, {"P_ox_tank": results["P_ox_tank"]}, validation_data, "Tank Pressure", "Pressure [Pa]")

            plt.subplot(1, 3, 2)
            plot_series(time, {"m_ox_tank": results["m_ox_tank"]}, validation_data, "Tank Mass", "Mass [kg]")

            plt.subplot(1, 3, 3)
            plot_series(time, {"m_dot_ox": results["m_dot_ox"]}, validation_data, "Mass Flow Rate", "kg/s")

        elif inputs.ox_tank_model == 2:  # non-equilibrium
            plt.subplot(2, 2, 1)
            plot_series(time, {"P_ox_tank": results["P_ox_tank"]}, validation_data, "Tank Pressure", "Pressure [Pa]")

            plt.subplot(2, 2, 2)
            plot_series(time, {"m_ox_tank": results["m_ox_tank"]}, validation_data, "Tank Mass", "Mass [kg]")

            plt.subplot(2, 2, 3)
            plot_series(time, {"m_dot_ox": results["m_dot_ox"]}, validation_data, "Mass Flow Rate", "kg/s")

            plt.subplot(2, 2, 4)
            plot_series(time, {"T_liq_ox_tank": results["T_liq_ox_tank"],
                               "T_sat_ox_tank": results["T_sat_ox_tank"],
                               "T_gas_ox_tank": results["T_gas_ox_tank"]
                               }, validation_data, "Temperature", "K")



    elif mode == "fuel_tank":
        plt.subplot(1, 1, 1)
        plot_series(time, {"P_fuel_tank": results["P_fuel_tank"]}, validation_data, "Fuel Tank Pressure", "Pressure [Pa]")

    elif mode == "full_stack":
        plt.subplot(1, 2, 1)
        plot_series(time, {"thrust": results["thrust"]}, validation_data, "Thrust Curve", "Thrust [N]")

        plt.subplot(1, 2, 2)
        plot_series(time, {
            "P_ox_tank": results["P_ox_tank"],
            "P_fuel_tank": results["P_fuel_tank"],
            "P_cc": results["P_cc"],
        }, validation_data, "System Pressures", "Pressure [Pa]")

    # Save/show
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()
