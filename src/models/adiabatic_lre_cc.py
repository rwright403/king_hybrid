###NOTE: THIS IS BAD! should have neglected mixing of helium and GOX!!!

from rocketcea.cea_obj_w_units import CEA_Obj #how to access this as a non us citizen?
import CoolProp.CoolProp as CP #I love coolprop! ~ units: http://www.coolprop.org/v4/apidoc/CoolProp.html
import matplotlib.pyplot as plt
import numpy as np

#TODO: DOUBLE CHECK PropsSI units       

class model():
    def __init__(self, oxidizer_name, fuel_name, A_throat, A_exit, P_atm, TIMESTEP):

        #use SI UNITS
        self.C = CEA_Obj(oxName=oxidizer_name, fuelName=fuel_name, pressure_units='Pa', isp_units='sec', cstar_units='m/s', temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='kJ/kg', density_units='kg/m^3', specific_heat_units='kJ/kg-K')
        #NOTE: SEEMS TO BREAK ABOVE

        self.P_cc = 0
        self.OF = 0

        self.y = 0 
        self.R = 0 # J/Kg*K

        # initial values (treat as private members)
        self.P_cc = P_atm
        self.P_atm = P_atm
        self.v_exit = 0
        self.r_dot_t = 0.1 

        self.m_dot_cc_t = 0        
        self.instThrust = 0

        self.A_throat = A_throat 
        self.A_exit = A_exit 
        self.expratio = self.A_exit / self.A_throat

        self.TIMESTEP = TIMESTEP

    #NOTE: no injector term currently in this script!!!!!
    def inst(self, m_dot_ox, m_dot_fuel):
        
        #print("this should not be nan: ",m_dot_ox,m_dot_fuel)
        #update O/F ratio and m_dot_cc
        self.m_dot_cc_t = m_dot_ox + m_dot_fuel
        self.OF = m_dot_ox/m_dot_fuel
        
        #CALL CEA to solve combustion
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K)
        self.y = fluid_prop[1] # (-)
        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.expratio, 0, 1)
        T_cc = temperatures[0]

        #print("TEMPERATURE: ",T_cc)
        
        #print(self.m_dot_cc_t, m_dot_fuel, m_dot_ox)

        #NOTE: P_cc too high?
        #print("AAAAA:", self.m_dot_cc_t, self.A_throat, self.R, T_cc, self.y)
        self.P_cc = (self.m_dot_cc_t / self.A_throat ) * np.sqrt( self.R*T_cc )  / ( np.sqrt( self.y * (2 / (self.y+1))**( (self.y+1)/(self.y-1) ) ) )

        #print("in cc:", self.P_cc)

        #resolve fluid properties
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K)
        self.y = fluid_prop[1] # (-)
        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.expratio, 0, 1)
        T_cc = temperatures[0]
        #T_throat = temperatures[1]
        #T_exit = temperatures[2]

        ###START: NOZZLE
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
        

        #solve thrust
        self.instThrust = (self.m_dot_cc_t * self.v_exit) + self.A_exit * (P_exit - self.P_atm)