#https://rocketcea.readthedocs.io/en/latest/functions.html#module-rocketcea.cea_obj

from numpy.random.mtrand import gamma
import numpy as np
import matplotlib.pyplot as plt
from src.models.cc._base import BaseChamber



class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, L_star, m_fuel_i, rho_fuel, a, n, L, A_port, nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        self.C = C
        #initalizing y and R
        self.y = 0 
        self.R = 0 # J/Kg*K

        self.rho_fuel = rho_fuel 
        self.m_fuel_t = m_fuel_i 
        self.a = a 
        self.n = n 

        #fuel grain geometry
        self.L = L
        self.A_port = A_port

        self.expratio = self.A_exit / self.A_throat
        #print(self.expratio)

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
        self.radius = np.sqrt(A_port / np.pi)

        # setup
        self.m_dot_fuel = 0
        self.total_propellant = 0

    def inst(self, m_dot_ox, _):
        # print( "ox m dot: ", m_dot_ox, " cc mass flow rate: ", self.m_dot_cc_t, " port area: ", self.A_port_t )
        # oxidizer mass flux iteration loop to get mass flow rate

        i = 0
        while i < 2:
            if i >= 1:
                G_ox_t = (m_dot_ox + self.m_dot_cc_t) / (2 * self.A_port)
            else:
                G_ox_t = m_dot_ox / self.A_port

            self.r_dot_t = self.a * (G_ox_t)**self.n
            
            #self.m_dot_fuel = self.rho_fuel * np.pi * ((self.r_dot_t + np.sqrt(self.A_port_t / np.pi))**2 - (self.A_port_t / np.pi)) * self.L
            #BUG: THE ABOVE IS WRONG, r_dot is bad!
            self.m_dot_fuel = self.rho_fuel * np.pi * ((self.r_dot_t*self.timestep + np.sqrt(self.A_port_t / np.pi))**2 - (self.A_port_t / np.pi)) * self.L
            
            if self.m_fuel_t <= 0:
                self.m_dot_fuel = 0
                self.OF = 0
            else:
                self.OF = m_dot_ox / self.m_dot_fuel

            i += 1

        self.m_dot_cc_t = m_dot_ox + self.m_dot_fuel

        # solve fluid properties
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K) #NOTE: IS THIS TOO BIG? UNIT ERROR RN?
        self.y = fluid_prop[1] # (-)

        #print(self.R, self.y)

        # call cea to get cc and exit temperature
        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.expratio, 0, 1)
        T_cc = temperatures[0]
        
        self.P_cc = (self.m_dot_cc_t /  self.A_throat ) * np.sqrt( self.R*T_cc )  / ( np.sqrt( self.y * (2 / (self.y+1))**( (self.y+1)/(self.y-1) ) ) )
        print(self.P_cc,self.m_dot_cc_t,self.A_throat,self.R,T_cc,self.y)


        #resolve fluid properties
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K)
        self.y = fluid_prop[1] # (-)

        # call cea to get cc and exit temperature
        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.expratio, 0, 1)
        T_cc = temperatures[0]

        #solve critical pressure:
        P_crit = self.P_cc*(1+ (self.y-1)/2)**((-1)*self.y/(self.y-1))


        
        #if subsonic at throat!
        if(P_crit < self.P_atm):
            #print("flow no longer choked! --> subsonic ")
            P_exit = self.P_atm
            T_exit = T_cc / (self.P_cc/P_exit)**((self.y-1)/self.y)
            #print(T_exit)
            cp = self.C.get_Chamber_Cp(self.P_cc, self.OF, self.expratio,0) #this is equilibrium check!
            self.v_exit = np.sqrt((2*cp)*(T_cc-T_exit))
            #print(self.v_exit, self.m_dot_cc_t)
        else:
        #supersonic #TODO: reverse case for speed? #TODO: make subsonic case a funciton where you pass in noteworthy upstream pressure
            if( 0.5 * self.P_atm < P_crit < 1.5* self.P_atm): #check that this is about equal and 50% deviation is about right
                #print("sonic throat subsonic exit case")
                P_exit = self.P_atm
                T_exit = T_cc / (P_crit/P_exit)**((self.y-1)/self.y) #pressure at throat will decrease from subsonic nozzle
                cp = self.C.get_Chamber_Cp(self.P_cc, self.OF, self.expratio,0) #this is equilibrium check!
                self.v_exit = np.sqrt((2*cp)*(T_cc-T_exit))

            else:
                exit_mach = self.C.get_MachNumber(self.P_cc, self.OF, self.expratio, 0, 1)

                P_exit = self.P_cc / (1 + ((self.y - 1) / 2) * (exit_mach**2) )**(self.y / (self.y - 1))

                self.v_exit = np.sqrt(((2 * self.y) / (self.y - 1)) * self.R * T_cc * (1 - (P_exit / self.P_cc)**((self.y - 1) / self.y)))
        



        
        self.instThrust = (self.m_dot_cc_t * self.v_exit) + self.A_exit * (P_exit - self.P_atm)


        # Recalculate new fuel grain size
        self.m_fuel_t = self.m_fuel_t - self.m_dot_fuel * self.timestep

        self.total_propellant = self.total_propellant + self.m_dot_cc_t * self.timestep

        self.radius = self.radius + self.r_dot_t * self.timestep
        self.A_port_t = np.pi * self.radius**2

        return {"P_cc": self.P_cc, "thrust": self.instThrust}



