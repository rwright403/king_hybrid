import numpy as np
from scipy.optimize import brentq
from src.models.nozzle._base import BaseNozzle

def area_mach_func(M, gamma, A_ratio):
    term1 = 1 / M
    term2 = (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
    term3 = (1 + (gamma-1)/2 * M**2)**((gamma+1)/(2*(gamma-1)))
    return (term1 * term2 * term3) - A_ratio


class basic_nozzle_model(BaseNozzle):
    def __init__(self, d_throat: float, expratio: float, P_atm: float, C):
        super().__init__(d_throat, expratio, P_atm, C)
    
    def sol_thrust(self, P_cc, T_cc, gamma, R):

        Critical_PR = ( 2/(gamma+1) ) ** ( -gamma/(gamma+1) )
        Actual_PR = self.P_atm/P_cc

        #properties at nozzle exit plane
        P_exit = 0.0
        v_exit = 0.0
        m_dot_exit = 0.0

        # Choking at throat
        if Actual_PR <= Critical_PR: # Find exit mass flow rate:
            M_exit = brentq(area_mach_func, 1.0001, 6, args=(gamma, self.A_exit/self.A_throat) ) # brent's method between Ma 1 and Ma 6             
            m_dot_exit = self.A_throat * P_cc * np.sqrt(gamma/(R*T_cc)) * ( (gamma+1)/2 )**(-(gamma+1)/(2*(gamma-1)))

        
        # Subsonic at throat
        else: # Find exit mass flow rate:
            M_exit = brentq(area_mach_func, 0.0001, .99999, args=(gamma, self.A_exit/self.A_throat) ) # brent's method between Ma 1 and Ma 6
            m_dot_exit = self.A_exit * P_cc * np.sqrt(gamma/(R*T_cc)) * (1 + (gamma-1)/2*M_exit**2)**(-(gamma+1)/(2*(gamma-1)))
        
        T_exit = T_cc / (1 +(gamma-1)/2*M_exit**2) # Approximating T_0 = T_cc

        a_exit = np.sqrt(R*gamma*T_exit)

        v_exit = a_exit*M_exit
        
        P_exit = P_cc / ((1 + (gamma-1)/2 * M_exit**2)**(gamma/(gamma-1)) )    

        self.instThrust = m_dot_exit*v_exit + self.A_exit*(P_exit - self.P_atm)

        return self.instThrust, m_dot_exit

"""
if self.P_atm < P_exit: # Underexpanded flow
if self.P_atm = P_exit: # Perfectly Expanded Flow
if self.P_atm > P_exit: # Overexpanded flow
"""