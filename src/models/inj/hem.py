import numpy as np
import CoolProp.CoolProp as CP
from src.models.inj._base import BaseInjector

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
            s_2 = CP.PropsSI('S', 'H', h_1, 'P', P_1, "N2O") #assuming isentropic, upstream entropy equals downstream entropy
            h_2_hem = CP.PropsSI('H', 'S', s_2, 'P', pres, "N2O")
            rho_2_hem = CP.PropsSI('D', 'S', s_2, 'P', pres, "N2O")
                
            m_dot_hem = self.C_inj * rho_2_hem * np.sqrt( 2 * np.abs(h_1 -  h_2_hem) )
                
            m_dot_hem_arr.append(m_dot_hem)

        m_dot_hem_crit = np.max(m_dot_hem_arr)
        P_crit = downstream_pres_arr[np.argmax(m_dot_hem_arr)]
        
        if P_2 < P_crit: #flow is choked
            m_dot_hem = m_dot_hem_crit
            return m_dot_hem

        else: #flow is unchoked

            s_1 = CP.PropsSI('S', 'H', h_1, 'P', P_1, "N2O")
            h_2_hem = CP.PropsSI('H', 'S', s_1, 'P', P_2, "N2O")
            rho_2_hem = CP.PropsSI('D', 'S', s_1, 'P', P_2, "N2O")
            m_dot_hem = self.C_inj * rho_2_hem * np.sqrt( 2 * (h_1 -  h_2_hem) )

            return m_dot_hem