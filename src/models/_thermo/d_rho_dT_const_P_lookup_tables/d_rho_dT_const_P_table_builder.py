import numpy as np
from thermo import Chemical
from src.models._thermo.n2o_thermo_span_wagner_class import SpanWagnerEOS_SingleState
import os
n2o = Chemical('nitrous oxide', T=298.15)
rho_crit = n2o.rhoc

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
gas_out_file = os.path.join(SCRIPT_DIR, "n2o_d_rho_dT_const_P_gas_table.npz")
liq_out_file = os.path.join(SCRIPT_DIR, "n2o_d_rho_dT_const_P_liq_table.npz")

def build_n2o_d_rho_dT_const_P_gas_table(
    rho_min=1, rho_max=rho_crit, nrho=200,
    T_min=180.0, T_max=330.0, nT=200,
    gas_out_file=gas_out_file
):
    rho_grid = np.linspace(rho_min, rho_max, nrho)   # Pa
    T_grid = np.linspace(T_min, T_max, nT)   # K

    d_rho_dT_const_P_gas_tab = np.empty((nrho, nT), dtype=np.float64)

    for i, rho in enumerate(rho_grid):
        for j, T in enumerate(T_grid):
            try:
                SWEOS = SpanWagnerEOS_SingleState(rho, T)
                d_rho_dT_const_P_gas_tab[i, j] = SWEOS.d_rho_dT_const_P   # W/m/K
            except Exception:
                d_rho_dT_const_P_gas_tab[i, j] = np.nan

    np.savez(gas_out_file, T_grid=T_grid, rho_grid=rho_grid, d_rho_dT_const_P_gas_tab=d_rho_dT_const_P_gas_tab)
    print("saved to:", gas_out_file)
    return rho_grid, T_grid, d_rho_dT_const_P_gas_tab

def build_n2o_d_rho_dT_const_P_liq_table(
    rho_min=rho_crit, rho_max=1300, nrho=200,
    T_min=180.0, T_max=330.0, nT=200,
    liq_out_file=liq_out_file
):
    rho_grid = np.linspace(rho_min, rho_max, nrho)   # Pa
    T_grid = np.linspace(T_min, T_max, nT)   # K

    d_rho_dT_const_P_liq_tab = np.empty((nrho, nT), dtype=np.float64)

    for i, rho in enumerate(rho_grid):
        for j, T in enumerate(T_grid):
            try:
                SWEOS = SpanWagnerEOS_SingleState(rho, T)
                d_rho_dT_const_P_liq_tab[i, j] = SWEOS.d_rho_dT_const_P   # W/m/K
            except Exception:
                d_rho_dT_const_P_liq_tab[i, j] = np.nan

    np.savez(liq_out_file, T_grid=T_grid, rho_grid=rho_grid, d_rho_dT_const_P_liq_tab=d_rho_dT_const_P_liq_tab)
    print("saved to:", liq_out_file)
    return rho_grid, T_grid, d_rho_dT_const_P_liq_tab

build_n2o_d_rho_dT_const_P_gas_table()
build_n2o_d_rho_dT_const_P_liq_table()