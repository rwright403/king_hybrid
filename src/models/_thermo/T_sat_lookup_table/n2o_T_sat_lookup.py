import numpy as np
import CoolProp.CoolProp as CP

class N2OSatTemperatureTable:
    def __init__(self, path="src/models/_thermo/T_sat_lookup_table/n2o_T_sat_table.npz"):
        d = np.load(path)
        self.P_grid = d["P_grid"]
        self.T_sat = d["T_sat"]
        self.logP_grid = np.log(self.P_grid)

    def lookup(self, P):
        if P < self.P_grid[0] or P > self.P_grid[-1]:
            # fallback if outside table
            print(f"Caution: lookup T_sat returned nan at (P={P:.0f} Pa), using lib directly (slow)")
            return CP.PropsSI("T", "P", P, "Q", 1, "N2O")

        return np.interp(np.log(P), self.logP_grid, self.T_sat)