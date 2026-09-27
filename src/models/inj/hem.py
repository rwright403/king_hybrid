import numpy as np
import CoolProp.CoolProp as CP
from src.models.inj._base import BaseInjector

from src.models._thermo.s_lookup_tables_f_h_P.n2o_s_lookup_f_h_P import N2O_s_f_h_P_Table
s_f_h_P_table = N2O_s_f_h_P_Table()
from src.models._thermo.h_lookup_tables_f_s_P.n2o_h_lookup_f_s_P import N2O_h_f_s_P_Table
h_f_s_P_table = N2O_h_f_s_P_Table()
from src.models._thermo.rho_lookup_tables_f_s_P.n2o_rho_lookup_f_s_P import N2O_rho_f_s_P_Table
rho_f_s_P_table = N2O_rho_f_s_P_Table()


class hem_model(BaseInjector):     
    """
    HEM MODEL with flow choking
    """
    def __init__(self, Cd: float, A_inj: float):
        super().__init__(Cd, A_inj)


    def m_dot(self, state: dict) -> float:
        
        P_1 = state["P_1"]
        P_2 = state["P_2"]
        h_1 = state["h_1"]

        downstream_pres_arr = np.linspace(P_2, P_1, 100)
        m_dot_hem_arr = []

        for pres in downstream_pres_arr:
            #s_2 = CP.PropsSI('S', 'H', h_1, 'P', P_1, "N2O") #assuming isentropic, upstream entropy equals downstream entropy
            s_2 = s_f_h_P_table.lookup(h_1, P_1)
            #h_2_hem = CP.PropsSI('H', 'S', s_2, 'P', pres, "N2O")
            h_2_hem = h_f_s_P_table.lookup(s_2, pres)
            #rho_2_hem = CP.PropsSI('D', 'S', s_2, 'P', pres, "N2O")
            rho_2_hem = rho_f_s_P_table.lookup(s_2, pres)

            m_dot_hem = self.C_inj * rho_2_hem * np.sqrt( 2 * np.abs(h_1 -  h_2_hem) )
                
            m_dot_hem_arr.append(m_dot_hem)

        m_dot_hem_crit = np.max(m_dot_hem_arr)
        P_crit = downstream_pres_arr[np.argmax(m_dot_hem_arr)]
        
        if P_2 < P_crit: #flow is choked
            m_dot_hem = m_dot_hem_crit
            return m_dot_hem

        else: #flow is unchoked

            #s_1 = CP.PropsSI('S', 'H', h_1, 'P', P_1, "N2O")
            s_1 = s_f_h_P_table.lookup(h_1, P_1)
            #h_2_hem = CP.PropsSI('H', 'S', s_1, 'P', P_2, "N2O")
            h_2_hem = h_f_s_P_table.lookup(s_1, P_2)
            #rho_2_hem = CP.PropsSI('D', 'S', s_1, 'P', P_2, "N2O")
            rho_2_hem = rho_f_s_P_table.lookup(s_1, P_2)

            m_dot_hem = self.C_inj * rho_2_hem * np.sqrt( 2 * (h_1 -  h_2_hem) )

            return m_dot_hem