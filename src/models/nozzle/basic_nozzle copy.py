import numpy as np
from src.models.nozzle._base import BaseNozzle

class basic_nozzle_model(BaseNozzle):
    def __init__(self, d_throat: float, expratio: float, P_atm: float, C):
        super().__init__(d_throat, expratio, P_atm, C)

    def sol_thrust(self, P_cc, T_cc, gamma, R):
        """
        Super-basic 0D nozzle model.
        Assumes always choked at throat and isentropic expansion to exit.
        """

        # --- 1. Mass flow rate (choked-flow formula)
        m_dot = (self.A_throat * P_cc / np.sqrt(T_cc)) \
              * np.sqrt(gamma / R) \
              * ( (2 / (gamma + 1)) ** ( (gamma + 1) / (2 * (gamma - 1)) ) )

        # --- 2. Exit conditions (isentropic expansion, always choked)
        P_exit = self.P_atm
        v_exit = np.sqrt( 2 * gamma / (gamma - 1) * R * T_cc *
                          (1 - (P_exit / P_cc) ** ((gamma - 1) / gamma)) )

        # --- 3. Thrust
        thrust = m_dot * v_exit + self.A_exit * (P_exit - self.P_atm)

        self.instThrust = thrust
        return thrust, m_dot
