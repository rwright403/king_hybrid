#https://rocketcea.readthedocs.io/en/latest/functions.html#module-rocketcea.cea_obj

from numpy.random.mtrand import gamma
import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
from src.models.cc._base import BaseChamber

R_UNIV = 8314 # J / (kg K)

class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, L_star, m_fuel_i, rho_fuel, a, n, L, A_port_i, nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        #use SI UNITS
        self.C = C

        self.OF = 0.5

        # fuel grain material
        self.rho_fuel = rho_fuel 
        self.a = a 
        self.n = n 

        #fuel grain geometry and mass
        self.m_fuel_t = m_fuel_i 
        self.L = L
        self.A_port_t = A_port_i
        self.radius = np.sqrt(A_port_i / np.pi)

        # initial values
        self.P_cc = P_atm
        self.P_atm = P_atm
        self.v_exit = 0
        self.r_dot_t = 0.1 

        self.m_cc = 0
        self.m_dot_reactants_out = 0   

        self.V_cc = L_star * self.nozzle.A_throat

        print("init check vol: ", self.V_cc, L_star, self.nozzle.A_throat)

        self.timestep = timestep    

        # setup
        self.m_dot_fuel = 0
        self.total_propellant = 0

    def inst(self, m_dot_ox, _: float = 0.0): #NOTE:m_dot_fuel solved from fuel grain in this cc. This is an artifact so program is modular
        # print( "ox m dot: ", m_dot_ox, " cc mass flow rate: ", self.m_dot_cc_t, " port area: ", self.A_port_t )
        # oxidizer mass flux iteration loop to get mass flow rate

        i = 0
        while i < 2:
            if i >= 1:
                G_ox_t = (m_dot_ox + self.m_dot_cc_t) / (2 * self.A_port_t)
            else:
                G_ox_t = m_dot_ox / self.A_port_t

            self.r_dot_t = self.a * (G_ox_t)**self.n
            
            self.m_dot_fuel = self.rho_fuel * np.pi * ((self.r_dot_t + np.sqrt(self.A_port_t / np.pi))**2 - (self.A_port_t / np.pi)) * self.L
            
            if self.m_fuel_t <= 0:
                self.m_dot_fuel = 0
                self.OF = 0
            else:
                self.OF = m_dot_ox / self.m_dot_fuel

            i += 1

            """fixed point iteration on P_cc"""

        for _ in range(10):

            #Update mass balance: 
            OF = m_dot_ox / self.m_dot_fuel
            m_dot_propellants_in = m_dot_ox + self.m_dot_fuel

            self.m_cc += (m_dot_propellants_in-self.m_dot_reactants_out)*self.timestep

            # CEA to solve combustion properties
            MW, y = self.C.get_Chamber_MolWt_gamma(self.P_cc, OF, self.nozzle.expratio)
            T_cc = self.C.get_Tcomb(self.P_cc, OF)

            # Treat reactant fluid as ideal gas and resolve P_cc
            R_spec = (R_UNIV/MW)
            self.P_cc = R_spec * (T_cc*self.m_cc/self.V_cc)

            #print("P_cc: ", self.P_cc, R_spec, T_cc, self.V_cc)

            instThrust, self.m_dot_reactants_out = self.nozzle.sol_thrust(self.P_cc, T_cc, y, R_spec)

            # As it loops it will update mass balance

        ### Recalculate new fuel grain size for next timestep
        self.m_fuel_t = self.m_fuel_t - self.m_dot_fuel * self.timestep

        self.radius = self.radius + self.r_dot_t * self.timestep
        self.A_port_t = np.pi * self.radius**2

        return {"P_cc": self.P_cc, "thrust": instThrust}

        """
        # solve fluid properties
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.nozzle.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K) #NOTE: IS THIS TOO BIG? UNIT ERROR RN?
        self.y = fluid_prop[1] # (-)

        #print(self.R, self.y)

        # call cea to get cc and exit temperature
        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.nozzle.expratio, 0, 1)
        T_cc = temperatures[0]
        
        # Use secant method to find chamber pressure
        #this only helps the first timestep otherwise useless
        #initial_enthalpy = self.C.get_Chamber_H(self.P_cc, self.OF, self.expratio)
        self.P_cc = (self.m_dot_cc_t / self.nozzle.A_throat ) * np.sqrt( self.R*T_cc ) / ( np.sqrt( self.y * (2 / (self.y+1))**( (self.y+1)/(self.y-1) ) ) )

        print("self.P_cc: ", self.P_cc, self.m_dot_cc_t , self.nozzle.A_throat, self.R, T_cc, self.y )        
        #print(self.P_cc,self.m_dot_cc_t,self.A_throat,self.R,T_cc,self.y)
        #while np.abs(merror(self.C, self.P_cc, self.OF, self.nozzle.expratio, self.nozzle.A_exit, self.m_dot_cc_t)) > 0.005:
        #   self.P_cc = secant(lambda P: merror(self.C, P, self.OF, self.nozzle.expratio, self.nozzle.A_exit, self.m_dot_cc_t), self.P_cc)
        #print(self.P_cc,self.m_dot_cc_t)
        #final_enthalpy = self.C.get_Chamber_H(self.P_cc, self.OF, self.expratio)

        #print(initial_enthalpy,final_enthalpy,initial_enthalpy-final_enthalpy)


        #resolve fluid properties
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.nozzle.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K)
        self.y = fluid_prop[1] # (-)

        # call CEA to get cc and exit temperature
        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.nozzle.expratio, 0, 1)
        T_cc = temperatures[0]
        #T_throat = temperatures[1]
        #T_exit = temperatures[2]
       
        self.instThrust = self.nozzle.sol_thrust( self.P_cc, T_cc, self.y, self.R)


        ### Recalculate new fuel grain size
        self.m_fuel_t = self.m_fuel_t - self.m_dot_fuel * self.timestep

        self.total_propellant = self.total_propellant + self.m_dot_cc_t * self.timestep

        self.radius = self.radius + self.r_dot_t * self.timestep
        self.A_port_t = np.pi * self.radius**2

        return {"P_cc": self.P_cc, "thrust": self.instThrust}




class model():
    def __init__(self, oxidizer_name, fuel_name, fuel_properties, m_fuel_i, rho_fuel, a, n, L, A_port_i, P_atm, A_throat, A_exit, timestep):
       
"""