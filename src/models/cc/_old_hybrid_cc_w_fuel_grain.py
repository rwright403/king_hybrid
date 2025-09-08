#https://rocketcea.readthedocs.io/en/latest/functions.html#module-rocketcea.cea_obj

from numpy.random.mtrand import gamma
import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
from src.models.cc._base import BaseChamber

from src.utils.math import secant

# Function to calculate the error in mass flow rate
def merror(C, P_cc, OF, expratio, A_exit, m_dot_in):
    #print("INSIDE MERROR", P_cc, OF, expratio, A_exit, m_dot_in)

    # Resolve CEA for given P_cc
    fluid_prop = C.get_Chamber_MolWt_gamma(P_cc, OF, expratio)
    R = 8314 / fluid_prop[0]  # kJ/(kg K)
    y = fluid_prop[1]  # (-)

    # Call CEA to get temperatures
    temperatures = C.get_Temperatures(P_cc, OF, expratio, 0, 1)
    T_exit = temperatures[2]

    # Call CEA to get densities
    densities = C.get_Densities(P_cc, OF, expratio, 0, 1)
    rho_exit = densities[2]

    # Solve for exit Mach number
    exit_mach = C.get_MachNumber(P_cc, OF, expratio, 0, 1)

    # Calculate exit velocity
    v_exit = exit_mach * np.sqrt(y * R * T_exit)

    # Calculate mass flow rate
    m_dot_i = rho_exit * v_exit * A_exit

    # Return the difference between calculated and target mass flow rate
    #print(m_dot_i, m_dot_in,m_dot_i - m_dot_in)
    return m_dot_i - m_dot_in



class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, m_fuel_i, rho_fuel, a, n, L, A_port_i, nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        #initalizing y and R
        self.y = 0 
        self.R = 0 # J/Kg*K

        self.rho_fuel = rho_fuel 
        self.m_fuel_t = m_fuel_i 
        self.a = a 
        self.n = n 

        #fuel grain geometry
        self.L = L
        self.A_port_t = A_port_i

        self.timestep = timestep

        # initial values (treat as private members)
        self.P_cc = P_atm
        self.P_atm = P_atm
        self.v_exit = 0
        self.r_dot_t = 0.1 

        self.m_floss = 0
        self.m_dot_cc_t = 0        
        self.instThrust = 0

        self.OF = 0.5
        self.radius = np.sqrt(A_port_i / np.pi)

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

        self.m_dot_cc_t = m_dot_ox + self.m_dot_fuel

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
       
        self.instThrust = self.nozzle.sol_thrust( P_cc=self.P_cc, T_cc=T_cc, m_dot=self.m_dot_cc_t, gamma=self.y, R=self.R)


        ### Recalculate new fuel grain size
        self.m_fuel_t = self.m_fuel_t - self.m_dot_fuel * self.timestep

        self.total_propellant = self.total_propellant + self.m_dot_cc_t * self.timestep

        self.radius = self.radius + self.r_dot_t * self.timestep
        self.A_port_t = np.pi * self.radius**2

        return {"P_cc": self.P_cc, "thrust": self.instThrust}





"""







class model():
    def __init__(self, oxidizer_name, fuel_name, fuel_properties, m_fuel_i, rho_fuel, a, n, L, A_port_i, P_atm, A_throat, A_exit, timestep):
       
"""