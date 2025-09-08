# src/models/tanks/_base.py
from abc import ABC, abstractmethod

class BaseTank(ABC):
    def __init__(self, injector, timestep: float):
        self.injector = injector
        self.timestep = timestep

    @abstractmethod
    def inst(self, P_cc: float):
        pass
