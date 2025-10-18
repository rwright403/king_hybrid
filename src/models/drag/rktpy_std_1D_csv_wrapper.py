import numpy as np
import os

def make_drag_func(path, mach_col=None, cd_col=None, delimiter=","):
    """
    Returns a function f(Mach, u, mu, rho) -> Cd that looks up/interpolates Cd vs Mach
    from a CSV file. CSV may have a header with 'Mach' and 'Cd' column names.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(path)

    # Read header to auto-detect columns if not provided
    with open(path, "r") as f:
        header = f.readline().strip()

    skip = 0
    if any(k in header.lower() for k in ("mach", "cd")):
        skip = 1
        names = [h.strip().lower() for h in header.split(delimiter)]
        if mach_col is None:
            # prefer exact 'mach' else first column
            mach_col = names.index("mach") if "mach" in names else 0
        if cd_col is None:
            cd_col = names.index("cd") if "cd" in names else len(names) - 1
    else:
        # no header: default to first and last columns
        if mach_col is None: mach_col = 0
        if cd_col is None:   cd_col   = -1

    data = np.loadtxt(path, delimiter=delimiter, skiprows=skip, ndmin=2)
    Mach = data[:, mach_col].astype(float)
    Cd   = data[:, cd_col].astype(float)

    # sort just in case
    idx = np.argsort(Mach)
    Mach, Cd = Mach[idx], Cd[idx]

    # build the callable RocketPy expects
    def cd_fn(mach, u, mu, rho):
        # 1-D linear interpolation in Mach; clamp at ends
        return float(np.interp(mach, Mach, Cd, left=Cd[0], right=Cd[-1]))

    return cd_fn