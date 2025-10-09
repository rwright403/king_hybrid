import os
import pandas as pd

def build_rocketpy_input_csv(df: pd.DataFrame, key: str, time_col="time", output_dir="src/tmp"):
    """
    Saves a CSV for a given key (e.g. 'm_dot_ox') from a DataFrame.
    Creates the output directory if it doesn't exist.
    """

    if key not in df.columns:
        raise ValueError(f"Key:'{key}', not found in DataFrame")

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{key}.csv")

    pd.DataFrame({
        "time": df[time_col],
        "value": df[key]
    }).to_csv(output_path, index=False)

    print(f"Saved {key} data to: {output_path}")
