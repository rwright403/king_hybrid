import numpy as np
import CoolProp.CoolProp as CP
from src.models.inj._base import BaseInjector

class spi_model(BaseInjector):
    def __init__(self, Cd: float, A_inj: float):
        super().__init__(Cd, A_inj)

    def m_dot(self, state: dict) -> float:
        P_1 = state["P_1"]
        P_2 = state["P_2"]
        rho_1 = state["rho_1"]
        
        m_dot_spi = self.C_inj * np.sqrt( 2 * rho_1 * (P_1-P_2) )

        return m_dot_spi