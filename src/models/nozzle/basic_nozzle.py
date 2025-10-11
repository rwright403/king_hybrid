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
        #NOTE: EQNS SET UP ACCORDING TO NOZZLE THEORY AND THERMODYNAMIC RELATIONS CH IN RPE
        # Assuming chamber properties are stagnation properties

        Critical_PR = ( (gamma+1)/2 ) ** ( gamma/(gamma-1) )
        Actual_PR = P_cc/self.P_atm

        #print("crit pr, actual pr: ", Critical_PR, Actual_PR)

        #initialize properties at nozzle exit plane
        P_exit = 0.0
        v_exit = 0.0
        m_dot_exit = 0.0
        M_exit=0.0

        #print("nozzle PR: ", Actual_PR, Critical_PR, P_cc)
        
        # Choking at throat
        if Actual_PR >= Critical_PR: # Find exit mass flow rate:
            M_exit = brentq(area_mach_func, 1.0001, 6, args=(gamma, self.expratio) ) # brent's method between Ma 1 and Ma 6             
            m_dot_exit = self.A_throat * P_cc * gamma * np.sqrt( (2/(gamma+1))**((gamma+1)/(gamma-1)) ) / np.sqrt(gamma*R*T_cc)

                    
            T_exit = T_cc / (1 +(gamma-1)/2*M_exit**2) 
            a_exit = np.sqrt(R*gamma*T_exit)
            v_exit = a_exit*M_exit
            P_exit = P_cc / ((1 + (gamma-1)/2 * M_exit**2)**(gamma/(gamma-1)) )  

            #print("choked", M_exit, m_dot_exit)
        
        # Subsonic at throat
        else: # Find exit mass flow rate:
            M_exit = brentq(area_mach_func, 0.0001, .99999, args=(gamma, self.A_exit/self.A_throat) ) # brent's method between Ma 1 and Ma 6
            
            #print("into T_exit: ", T_cc, gamma, M_exit)

            T_exit = T_cc / (1 +(gamma-1)/2*M_exit**2)
            P_exit = P_cc / ((1 + (gamma-1)/2 * M_exit**2)**(gamma/(gamma-1))) 
            
            #print("denom: ", R, T_exit)
            
            rho_exit = P_exit / (R*T_exit)
            a_exit = np.sqrt(gamma*R*T_exit)
            v_exit = M_exit * a_exit

            m_dot_exit = rho_exit * v_exit * self.A_exit
            
            print("NOZZLE IS UNCHOKED AT THROAT", M_exit, m_dot_exit, T_exit)
 

        self.instThrust = m_dot_exit*v_exit + self.A_exit*(P_exit - self.P_atm)

        #print("nozzle: ", self.instThrust, m_dot_exit, v_exit, P_exit, self.P_atm) 

        return self.instThrust, m_dot_exit

        

"""
if self.P_atm < P_exit: # Underexpanded flow
if self.P_atm = P_exit: # Perfectly Expanded Flow
if self.P_atm > P_exit: # Overexpanded flow
"""