import os
import numpy as np
from src.models._thermo.n2o_thermo_span_wagner_class import SpanWagnerEOS_EquilibriumPhase

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def build_n2o_sat_T_table(
    P_min=9.0e4,
    P_max=7.2e6,
    nP=1000,
    out_file=None,
):
    if out_file is None:
        out_file = os.path.join(SCRIPT_DIR, "n2o_T_sat_table.npz")

    P_grid = np.geomspace(P_min, P_max, nP)   # Pa, log-spaced
    T_sat = np.empty_like(P_grid)

    for i, P in enumerate(P_grid):
        try:
            SWEOS_eq = SpanWagnerEOS_EquilibriumPhase(None, P)
            T_sat[i] = SWEOS_eq.T
        except Exception:
           T_sat[i] = np.nan

    np.savez(out_file, P_grid=P_grid, T_sat=T_sat)
    print("saved to:", out_file)
    return P_grid, T_sat

#build_n2o_sat_T_table()