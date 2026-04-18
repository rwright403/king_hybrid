import numpy as np
import CoolProp.CoolProp as CP

class N2OSatPressureTable:
    def __init__(self, path="src/models/_thermo/P_sat_lookup_table/n2o_P_sat_table.npz"):
        d = np.load(path)
        self.T_grid = d["T_grid"]
        self.P_sat = d["P_sat"]
        self.logT_grid = np.log(self.T_grid)

    def lookup(self, T):
        if T < self.T_grid[0] or T > self.T_grid[-1]:
            # fallback if outside table
            print(f"Caution: lookup T_sat returned nan at (T={T:.5f} K), using lib directly (slow)")
            return CP.PropsSI("P", "T", T, "Q", 1, "N2O")

        return np.interp(np.log(T), self.logT_grid, self.P_sat)