import numpy as np
import CoolProp.CoolProp as CP

class N2ORhoSatGasTable:
    def __init__(self, path="src/models/_thermo/rho_sat_lookup_tables/n2o_rho_sat_gas_table.npz"):
        d = np.load(path)
        self.P_grid = d["P_grid"]
        self.rho_sat_gas = d["rho_sat_gas"]
        self.logP_grid = np.log(self.P_grid)

    def lookup(self, P):
        return np.interp(np.log(P), self.logP_grid, self.rho_sat_gas)