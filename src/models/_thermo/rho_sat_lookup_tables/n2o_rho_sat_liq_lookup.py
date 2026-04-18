import numpy as np
import CoolProp.CoolProp as CP

class N2ORhoSatLiqTable:
    def __init__(self, path="src/models/_thermo/rho_sat_lookup_tables/n2o_rho_sat_liq_table.npz"):
        d = np.load(path)
        self.P_grid = d["P_grid"]
        self.rho_sat_liq = d["rho_sat_liq"]
        self.logP_grid = np.log(self.P_grid)

    def lookup(self, P):
        return np.interp(np.log(P), self.logP_grid, self.rho_sat_liq)