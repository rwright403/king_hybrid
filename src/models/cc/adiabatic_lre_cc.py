import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
from src.models.cc._base import BaseChamber

R_UNIV = 8314 # J / (kg K)

#TODO: DOUBLE CHECK PropsSI units  

class adiabatic_lre_cc_model(BaseChamber):

    def __init__(self, L_star, nozzle, P_atm, C, timestep):
        self.nozzle = nozzle

        #use SI UNITS
        self.C = C

        self.OF = 0.5

        # initial values
        self.P_atm = P_atm
        self.P_cc = P_atm
        self.v_exit = 0
        self.r_dot_t = 0.1 

        self.m_cc = 0
        self.m_dot_reactants_out = 0

        self.V_cc = L_star * self.nozzle.A_throat

        self.timestep = timestep


    def inst(self, m_dot_ox, m_dot_fuel):

        """fixed point iteration on P_cc"""

        for _ in range(10):

            #Update mass balance: 
            OF = m_dot_ox / m_dot_fuel
            m_dot_propellants_in = m_dot_ox + m_dot_fuel

            self.m_cc += (m_dot_propellants_in-self.m_dot_reactants_out)*self.timestep

            # CEA to solve combustion properties
            MW, y = self.C.get_Chamber_MolWt_gamma(self.P_cc, OF, self.nozzle.expratio)
            T_cc = self.C.get_Tcomb(self.P_cc, OF)

            # Treat reactant fluid as ideal gas and resolve P_cc
            R_spec = (R_UNIV/MW)
            self.P_cc = R_spec * (T_cc*self.m_cc/self.V_cc)

            instThrust, self.m_dot_reactants_out = self.nozzle.sol_thrust(self.P_cc, T_cc, y, R_spec)

            # As it loops it will update mass balance and converge
            #print(self.P_cc, instThrust)

        return {"P_cc": self.P_cc, "thrust": instThrust}