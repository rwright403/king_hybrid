#https://rocketcea.readthedocs.io/en/latest/functions.html#module-rocketcea.cea_obj

from numpy.random.mtrand import gamma
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel



class cc():
    def __init__(self, oxidizer_name, fuel_name, fuel_properties, m_fuel_i, rho_fuel, a, n, L,A_port_i, P_atm, A_throat, A_exit, timestep):
       
        ### Setup RocketCEA fuel characteristics

        add_new_fuel(fuel_name, fuel_properties)

        self.C = CEA_Obj(oxName=oxidizer_name, fuelName=fuel_name, pressure_units='Pa', isp_units='sec', cstar_units='m/s', temperature_units='K', sonic_velocity_units='m/s',enthalpy_units='kJ/kg',density_units='kg/m^3',specific_heat_units='kJ/kg-K')

        #fuel grain
        self.y = 1.23
        self.R = 300 #J/Kg*K

        self.rho_fuel = rho_fuel 
        self.m_fuel_t = m_fuel_i 
        self.a = a 
        self.n = n 

        #geometry
        self.L = L
        self.A_port_t = A_port_i

        self.A_throat = A_throat 
        self.A_exit = A_exit 

        self.timestep = timestep

        
        #intial values
        self.P_cc = P_atm
        ### UNUSED #self.A_port_f = 0.00411180831 #m^3
        self.r_dot_t = 0.1 




        self.m_floss = 0
        self.m_dot_cc_t = 0        
        self.instThrust = 0

        self.OF = 0.5
        self.radius = np.sqrt(A_port_i/np.pi)

        #setup
        self.m_dot_fuel = 0#self.L*self.rho_fuel*self.r_dot_t*(2*np.sqrt(self.A_port_t*np.pi))




        

    def inst(self, m_dot_ox, P_exit):


        #print( "ox m dot: ", m_dot_ox, " cc mass flow rate: ", self.m_dot_cc_t, " port area: ", self.A_port_t )
        #oxidizer mass flux iteration loop to get mass flow rate
        i = 0
        while i < 5:
            if i>1:
                G_ox_t = (m_dot_ox + self.m_dot_cc_t )/(2*self.A_port_t)
            else:
                G_ox_t = m_dot_ox/self.A_port_t
            i+=1

            #print("Gox: ", G_ox_t,  " m_dot_ox: ", m_dot_ox, "self.m_dot_cc_t", self.m_dot_cc_t, " self.A_port_t: ", self.A_port_t)

            self.r_dot_t = self.a*(G_ox_t)**self.n
            
            m_dot_fuel = self.L*self.rho_fuel*self.r_dot_t*(2*np.sqrt(self.A_port_t*np.pi))

            self.m_dot_cc_t  = m_dot_ox + m_dot_fuel

            i = i + 1

        

        #solve flame temperature
        T_cc = self.C.get_Tcomb(self.P_cc,(self.OF))
        #print("flame temp k: ", T_cc)

        #solve fluid properties
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, (self.A_exit/self.A_throat))
        self.R = 8314/fluid_prop[0]
        self.y = fluid_prop[1]

                       
        #solve exit velocity
        v_exit = np.sqrt( ((2*self.y)/(self.y-1)) *self.R*T_cc* (1-(P_exit/self.P_cc)**((self.y-1)/self.y) ) )
        #print(v_exit)

        #solve exit thrust
        self.instThrust = ((m_dot_fuel+m_dot_ox)*v_exit) #+ self.A_exit*(70000-101325)   #TODO: COMBINE W FLIGHT SIM MODEL

        #print("ISP: ", self.instThrust/(self.m_dot_cc_t*9.81))

        cstar = self.C.get_Cstar(self.P_cc, self.OF)
        #print("cstar: ", cstar)

        #recalculate P_cc
        self.P_cc = (self.m_dot_cc_t*cstar)/(2*self.A_throat) #this doesnt maek sense?
        #print("Pcc: ", self.P_cc, " m_dot_cc ", self.m_dot_cc_t, " cstar: ", cstar)

        #recalculate new fuel grain size
        self.m_floss = m_dot_fuel*self.timestep + self.m_floss

        self.radius = self.radius + self.r_dot_t * self.timestep

        self.A_port_t = np.pi *(self.radius**2)


        #recalculate O/F ratio
        self.OF = m_dot_ox/m_dot_fuel

        #print("Pcc: ", self.P_cc, " A_port_t: ", self.A_port_t, " m_dot_fuel: ", m_dot_fuel)