import numpy as np
from thermo import Chemical
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
gas_out_file = os.path.join(SCRIPT_DIR, "n2o_kg_table.npz")
liq_out_file = os.path.join(SCRIPT_DIR, "n2o_kl_table.npz")

def build_n2o_kg_table(
    T_min=180.0, T_max=330.0, nT=200,
    P_min=1.0e5, P_max=8.0e6, nP=200,
    gas_out_file=gas_out_file
):
    T_grid = np.linspace(T_min, T_max, nT)   # K
    P_grid = np.linspace(P_min, P_max, nP)   # Pa

    kg_tab = np.empty((nT, nP), dtype=np.float64)

    for i, T in enumerate(T_grid):
        for j, P in enumerate(P_grid):
            try:
                n2o = Chemical("N2O", T=T, P=P)
                kg_tab[i, j] = n2o.kg   # W/m/K
            except Exception:
                kg_tab[i, j] = np.nan

    np.savez(gas_out_file, T_grid=T_grid, P_grid=P_grid, kg_tab=kg_tab)
    print("saved to:", gas_out_file)
    return T_grid, P_grid, kg_tab

def build_n2o_kl_table(
    T_min=180.0, T_max=330.0, nT=200,
    P_min=1.0e5, P_max=8.0e6, nP=200,
    liq_out_file=liq_out_file
):
    T_grid = np.linspace(T_min, T_max, nT)   # K
    P_grid = np.linspace(P_min, P_max, nP)   # Pa

    kl_tab = np.empty((nT, nP), dtype=np.float64)

    for i, T in enumerate(T_grid):
        for j, P in enumerate(P_grid):
            try:
                n2o = Chemical("N2O", T=T, P=P)
                kl_tab[i, j] = n2o.kl   # W/m/K
            except Exception:
                kl_tab[i, j] = np.nan

    np.savez(liq_out_file, T_grid=T_grid, P_grid=P_grid, kl_tab=kl_tab)
    print("saved to:", liq_out_file)
    return T_grid, P_grid, kl_tab


#build_n2o_kl_table()
#build_n2o_kg_table()
