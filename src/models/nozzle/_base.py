import numpy as np
from abc import ABC, abstractmethod

class BaseNozzle(ABC):
    def __init__(self, d_throat: float, expratio: float, P_atm: float, C):
        self.C = C
        
        self.expratio = expratio
        self.A_throat = 0.25* np.pi * (d_throat)**2
        self.A_exit = expratio * self.A_throat

        self.P_atm = P_atm

        self.instThrust = 0

    @abstractmethod
    def sol_thrust(self, P_cc: float, T_cc: float, m_dot: float, gamma: float, R: float, chamber_obj=None):
        pass