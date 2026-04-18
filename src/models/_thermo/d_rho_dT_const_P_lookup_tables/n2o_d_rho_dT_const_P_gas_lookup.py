import numpy as np
from src.models._thermo.n2o_thermo_span_wagner_class import SpanWagnerEOS_SingleState

class N2ODRHODTPGasTable:
    def __init__(self, path="src/models/_thermo/d_rho_dT_const_P_lookup_tables/n2o_d_rho_dT_const_P_gas_table.npz"):
        d = np.load(path)
        self.rho_grid = d["rho_grid"]
        self.T_grid = d["T_grid"]
        self.d_rho_dT_const_P_gas_tab = d["d_rho_dT_const_P_gas_tab"]

    def lookup(self, rho, T):
        rho_grid = self.rho_grid
        T_grid = self.T_grid
        F = self.d_rho_dT_const_P_gas_tab

        i = np.searchsorted(rho_grid, rho) - 1 
        j = np.searchsorted(T_grid, T) - 1

        i = max(0, min(i, len(rho_grid) - 2))
        j = max(0, min(j, len(T_grid) - 2))        
   
        RHO1, RHO2 = rho_grid[i], rho_grid[i+1]
        T1, T2 = T_grid[j], T_grid[j+1]

        Q11 = F[i,   j]
        Q21 = F[i+1, j]
        Q12 = F[i,   j+1]
        Q22 = F[i+1, j+1]        

        # fallback if any NaN
        if np.isnan(Q11) or np.isnan(Q21) or np.isnan(Q12) or np.isnan(Q22):
            print(f"Caution: d_rho_dT_const_P_gas lookup returned nan at (rho={rho:.0f} kg/m^3), (T={T:.4f} K), using lib directly (approx .012s per call)")
            return SpanWagnerEOS_SingleState(rho, T).d_rho_dT_const_P
        
        a = (rho - RHO1)/(RHO2 - RHO1)
        b = (T - T1)/(T2 - T1)


        return (
            (1-a)*(1-b)*Q11 +
            a*(1-b)*Q21 +
            (1-a)*b*Q12 +
            a*b*Q22
        )