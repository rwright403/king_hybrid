from abc import ABC, abstractmethod

class BaseChamber(ABC):
    def __init__(self, nozzle, P_atm, timestep, C):
        self.nozzle = nozzle
        self.P_cc = P_atm
        self.P_atm = P_atm
        self.timestep = timestep
        self.C = C

        self.instThrust = 0.0
        self.m_dot_cc_t = 0.0

    @abstractmethod
    def inst(self, m_dot_ox, m_dot_fuel=0.0):
        pass
