from abc import ABC, abstractmethod

def build_state():
    state = {
        # Injector inlet (tank exit node)
        "P_1": None,              # Tank pressure at injector inlet [Pa]
        "P_2": None,              # Downstream (chamber) pressure [Pa]
        "T_1": None,              # Temperature at injector inlet [K]
        "rho_1": None,    # Density at injector inlet [kg/m^3]
        "h_1": None,              # Specific enthalpy at injector inlet [J/kg]
        "u_1": None,              # Specific internal energy at injector inlet [J/kg]
        "s_1": None,              # Specific entropy at injector inlet [J/kg-K]

        # Phase information
        "x_1": None,              # Vapor quality at injector inlet [-] (0 = liq, 1 = vap)
        "P_sat": None,            # Saturation pressure at tank T [Pa]
    }
    return state



class BaseInjector(ABC):
    def __init__(self, Cd: float, A_inj: float):
        self.C_inj = Cd*A_inj

    @abstractmethod
    def m_dot(self, state):
        """Return mass flow (kg/s) given upstream/downstream pressures and tank state."""
        pass
