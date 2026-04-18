import numpy as np
from thermo import Chemical

class N2OLiqConductivityTable:
    def __init__(self, path="src/models/_thermo/conductivity_lookup_tables/n2o_kl_table.npz"):
        d = np.load(path)
        self.T_grid = d["T_grid"]
        self.P_grid = d["P_grid"]
        self.kl_tab = d["kl_tab"]

    def lookup(self, T, P):
        T_grid = self.T_grid
        P_grid = self.P_grid
        F = self.kl_tab

        i = np.searchsorted(T_grid, T) - 1
        j = np.searchsorted(P_grid, P) - 1

        i = max(0, min(i, len(T_grid) - 2))
        j = max(0, min(j, len(P_grid) - 2))

        T1, T2 = T_grid[i], T_grid[i+1]
        P1, P2 = P_grid[j], P_grid[j+1]

        Q11 = F[i,   j]
        Q21 = F[i+1, j]
        Q12 = F[i,   j+1]
        Q22 = F[i+1, j+1]

        # fallback if any NaN
        if np.isnan(Q11) or np.isnan(Q21) or np.isnan(Q12) or np.isnan(Q22):
            print("Caution: kg lookup returned nan at (T={T:.4f} K), (P={P:.0f} Pa), using lib directly (approx .012s per call)")
            return Chemical("N2O", T=T, P=P).kl

        a = (T - T1)/(T2 - T1)
        b = (P - P1)/(P2 - P1)

        return (
            (1-a)*(1-b)*Q11 +
            a*(1-b)*Q21 +
            (1-a)*b*Q12 +
            a*b*Q22
        )