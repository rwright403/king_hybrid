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
        self.expratio = self.A_exit/self.A_throat

        self.timestep = timestep

        
        #intial values (treat as private members)
        self.P_cc = P_atm
        self.P_atm = P_atm
        self.v_exit = 0
        ### UNUSED #self.A_port_f = 0.00411180831 #m^3
        self.r_dot_t = 0.1 




        self.m_floss = 0
        self.m_dot_cc_t = 0        
        self.instThrust = 0

        self.OF = 0.5
        self.radius = np.sqrt(A_port_i/np.pi)

        #setup
        self.m_dot_fuel = 0#self.L*self.rho_fuel*self.r_dot_t*(2*np.sqrt(self.A_port_t*np.pi))




        

    def inst(self, m_dot_ox, COMBUSTION_PRESSURE):


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
            
            self.m_dot_fuel = self.rho_fuel*np.pi*( (self.r_dot_t+np.sqrt(self.A_port_t/np.pi) )**2 - (self.A_port_t/np.pi) )*self.L#self.L*self.rho_fuel*self.r_dot_t*(2*np.sqrt(self.A_port_t*np.pi))
            if(self.m_fuel_t<=0):
                self.m_dot_fuel = 0
                self.OF = 0

            self.m_dot_cc_t  = m_dot_ox + self.m_dot_fuel

            i += 1

        

        #solve flame temperature
        T_cc = self.C.get_Tcomb(self.P_cc,(self.OF))
        #print("flame temp k: ", T_cc, self.instThrust)

        #solve fluid properties
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = 8314/fluid_prop[0] #problem, this is likely dropping once either P_cc falls or OF falls!!!!!
        self.y = fluid_prop[1]

        #print(self.m_fuel_t, self.OF)

        #Solve Exit pressure
        exit_mach = self.C.get_MachNumber(self.P_cc, self.OF, self.expratio,1,0)
        #print(exit_mach, self.P_cc)

        P_exit = self.P_cc*(1+((self.y-1)/2)*exit_mach**2)**(self.y/(1-self.y))
        
        M_cc = self.C.get_Chamber_MachNumber(self.P_cc,self.OF,3) #TODO: ADD CR TO CONSTANTS!!!!
        T_0 = T_cc*(1+(self.y-1)/2*M_cc**2)

        #if we recalculate chamber pressure with cstar we get a more accurate chamber pressure.
        P_0= self.m_dot_cc_t/(self.A_throat*np.sqrt(self.y/(self.R*T_0))*(2/(self.y+1))**((self.y+1)/2*(self.y-1)))
        self.P_cc = P_0*(1+(self.y-1)/2*M_cc**2)**((-self.y)/(self.y-1))

        #cstar = self.C.get_Cstar(self.P_cc, self.OF)
        #self.P_cc=(self.m_dot_cc_t*cstar)/(self.A_throat)
                       
        #solve exit velocity
        self.v_exit = np.sqrt( ((2*self.y)/(self.y-1)) *self.R*T_cc* (1-(P_exit/self.P_cc)**((self.y-1)/self.y) ) )
        #print(self.v_exit," ",self.y," ",self.R," ",T_cc," ",P_exit," ",self.P_cc)

        #solve exit thrust 
        self.instThrust = (self.m_dot_cc_t*self.v_exit) +self.A_exit*(P_exit - self.P_atm)
        #print(self.v_exit,self.m_dot_cc_t, )
        #print("ISP: ", self.instThrust/(self.m_dot_cc_t*9.81))
        #print( self.m_dot_cc_t,self.v_exit,self.A_exit,P_exit,self.P_atm)

        if np.isnan(self.instThrust):
            self.instThrust = 0

        if self.instThrust == 0:
            self.instThrust = 0

        #if we recalculate chamber pressure with cstar we get a more accurate chamber pressure and a less accurate thrust
        cstar = self.C.get_Cstar(self.P_cc, self.OF)
        self.P_cc=(self.m_dot_cc_t*cstar)/(self.A_throat)

        #print("cstar: ", cstar, "P_cc: ", self.P_cc)


        print(P_0, self.P_cc, self.m_fuel_t)
    
            
        #recalculate new fuel grain size
        self.m_fuel_t = self.m_fuel_t-self.m_dot_fuel*self.timestep

        self.radius = self.radius + self.r_dot_t * self.timestep

        self.A_port_t = np.pi *(self.radius**2)


        #recalculate O/F ratio
        self.OF = m_dot_ox/self.m_dot_fuel




"""
%                                             &@&.                       
%                                          @@      /@                    
%                               %@@@@@@,  @&    @%   %(                  
%                           (@%         @@@        @                     
%              ,&@@@@@@@@@@@.         @@&         @#                     
%          *@@@@@@&      @/         @@,       ,&,  /@@@.                 
%         @@@@@%        @    &@@@@@@.                 @@%                
%        #@@@@@        @..@*    @@                     @@                
%        *@@@@@        @,    (@/                      &@,                
%         @@@@@@          @@.         *@@@@@,        #@#                 
%          @@@@@@    (@#           #@@      @       @@.                  
%            @@@@@@  .&@@@@@@    @@ @      @/     /@&                    
%             #@@@@@@.    #@   &@  @      @     @@/  #@,                 
%               .@@@@@@@. @@  @@@  @    @.   @@%     @@@%                
%               @  @@@@@@@@@ % @  ,   @%@@@*         #@@@                
%             /#      %@@@@@@@@@.                    @@@@/                       
%            /%           @@@@@@@@@@@@,           (@@@@@@                
%             @          *@.  *@@@@@@@@@@@@@@@@@@@@@@@@@                 
%            @/      .@@            ,&@@@@@@@@@@@@@@@                    
%           @    @@,                                                     
%          @@%                               
HAWK EM MAGPI!!!!"""