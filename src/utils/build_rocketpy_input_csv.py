import os
import pandas as pd

def build_rocketpy_input_csv(df: pd.DataFrame, key: str, time_col="time", output_dir="src/tmp"):
    """
    Saves a single [time, value] CSV for RocketPy from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing all simulation results.
    key : str
        The column name to export as 'value' (e.g., 'thrust', 'm_dot_ox').
    time_col : str
        The column representing time (default: 'time').
    output_dir : str
        Directory to save CSV (default: 'src/tmp').

    Returns
    -------
    str
        Path to the saved CSV file.
    """
    if time_col not in df.columns:
        raise ValueError(f"Time column '{time_col}' not found in DataFrame")
    if key not in df.columns:
        raise ValueError(f"Key:'{key}', not found in DataFrame")

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{key}.csv")

    print("output_path: ", output_path)

    pd.DataFrame({
        "time": df[time_col],
        "value": df[key]
    }).to_csv(output_path, index=False)

    return output_path

