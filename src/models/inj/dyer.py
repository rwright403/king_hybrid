import numpy as np
from src.models.inj._base import BaseInjector
from src.models.inj.spi import spi_model
from src.models.inj.hem import hem_model

class dyer_model(BaseInjector):

    def __init__(self, Cd: float, A_inj: float):
        super().__init__(Cd, A_inj)
        self.spi = spi_model(Cd, A_inj)
        self.hem = hem_model(Cd, A_inj)


    def m_dot(self, state: dict) -> float:

        P_1 = state["P_1"]
        P_2 = state["P_2"]
        P_sat = state["P_sat"]

        # First compute SPI
        m_dot_spi = self.spi.m_dot(state)
                        
        #NOTE: FOR THIS CASE NO CAVITATION SO WE ARE JUST USING THE SPI MODEL
        if(P_sat < P_2):
            m_dot_dyer = m_dot_spi

        #NOTE: ELSE TWO PHASE AT INJ OUTLET AND USE DYER TO ACCOUNT FOR TWO PHASE EFFECTS
        else:
            m_dot_hem = self.hem.m_dot(state)
            
            dyer_k = np.sqrt( (P_1 - P_2) / (P_sat - P_2) ) 
            m_dot_dyer = ((dyer_k/(1+dyer_k)) * m_dot_spi) + ((1/(1+dyer_k)) * m_dot_hem)
        
        return m_dot_dyer