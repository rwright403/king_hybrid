import os
import numpy as np
import CoolProp.CoolProp as CP

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def build_n2o_sat_T_table(
    T_min=180.23,
    T_max=309.52,
    nP=1000,
    out_file=None,
):
    if out_file is None:
        out_file = os.path.join(SCRIPT_DIR, "n2o_P_sat_table.npz")

    T_grid = np.geomspace(T_min, T_max, nP)   # Pa, log-spaced
    P_sat = np.empty_like(T_grid)

    for i, T in enumerate(T_grid):
        P_sat[i] = CP.PropsSI("P", "T", T, "Q", 1, "N2O")

    np.savez(out_file, T_grid=T_grid, P_sat=P_sat)
    print("saved to:", out_file)
    return T_grid, P_sat

#build_n2o_sat_T_table()