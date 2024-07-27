#https://rocketcea.readthedocs.io/en/latest/functions.html#module-rocketcea.cea_obj

from numpy.random.mtrand import gamma
import numpy as np
import matplotlib.pyplot as plt
#import CoolProp.CoolProp as CP
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel

# Define the secant method
def secant(func, x1, x_eps=0.005, kk_max=1000):
    x_eps = x1 * x_eps  # Set the tolerance to be a percentage of the initial guess
    x2 = x1 - x1 * 0.01  # Set a second point 1% away from the original guess
    F1 = func(x1)  # Evaluate function at x1
    F2 = func(x2)  # Evaluate function at x2
    #print("SECANT SETUP:", x1, x2, F1, F2)
    kk = 1  # Set up counter

    while np.abs(x2 - x1) >= x_eps and kk < kk_max:  # While error is too large and counter is less than max
        x3 = x2 - (F2 * (x2 - x1) / (F2 - F1))  # Secant method update
        x1 = x2  # Move everything forward
        x2 = x3
        F1 = F2
        F2 = func(x2)
        kk += 1
    return x2

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

class cc():
    def __init__(self, oxidizer_name, fuel_name, fuel_properties, m_fuel_i, rho_fuel, a, n, L, A_port_i, P_atm, A_throat, A_exit, timestep):
       
        ### Setup RocketCEA fuel characteristics

        add_new_fuel(fuel_name, fuel_properties)

        self.C = CEA_Obj(oxName=oxidizer_name, fuelName=fuel_name, pressure_units='Pa', isp_units='sec', cstar_units='m/s', temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='kJ/kg', density_units='kg/m^3', specific_heat_units='kJ/kg-K')

        # fuel grain
        #initalizing y and R
        self.y = 0 
        self.R = 0 # J/Kg*K

        self.rho_fuel = rho_fuel 
        self.m_fuel_t = m_fuel_i 
        self.a = a 
        self.n = n 

        # geometry
        self.L = L
        self.A_port_t = A_port_i

        self.A_throat = A_throat 
        self.A_exit = A_exit 
        self.expratio = self.A_exit / self.A_throat
        print(self.expratio)

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

    def inst(self, m_dot_ox):
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
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K)
        self.y = fluid_prop[1] # (-)

        #print(self.R, self.y)

        # call cea to get cc and exit temperature
        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.expratio, 0, 1)
        T_cc = temperatures[0]
        
        # Use secant method to find chamber pressure
        #this only helps the first timestep fucking useless
        #initial_enthalpy = self.C.get_Chamber_H(self.P_cc, self.OF, self.expratio)
        self.P_cc = (self.m_dot_cc_t /  self.A_throat ) * np.sqrt( self.R*T_cc )  / ( np.sqrt( self.y * (2 / (self.y+1))**( (self.y+1)/(self.y-1) ) ) )
        #while np.abs(merror(self.C, self.P_cc, self.OF, self.expratio, self.A_exit, self.m_dot_cc_t)) > 0.005:
         #   self.P_cc = secant(lambda P: merror(self.C, P, self.OF, self.expratio, self.A_exit, self.m_dot_cc_t), self.P_cc)
        #print(self.P_cc,self.m_dot_cc_t)
        #final_enthalpy = self.C.get_Chamber_H(self.P_cc, self.OF, self.expratio)

        #print(initial_enthalpy,final_enthalpy,initial_enthalpy-final_enthalpy)


        # solve fluid properties
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K)
        self.y = fluid_prop[1] # (-)

        # call cea to get cc and exit temperature
        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.expratio, 0, 1)
        T_cc = temperatures[0]
        #T_throat = temperatures[1]
        #T_exit = temperatures[2]

        #solve critical pressure:
        P_crit = self.P_cc*(1+ (self.y-1)/2)**((-1)*self.y/(self.y-1))

        #print(exit_mach,P_exit, P_crit, self.P_atm)

        #print(P_exit, "|", P_crit, self.P_atm,"|", exit_mach)
        
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
        
        #v_exit = exit_mach * np.sqrt(self.y*self.R*T_exit)


        #print(self.v_exit, v_exit)

        #print(P_exit)

        #this is wrong use v_cc

        #enthalpies =self.C.get_Enthalpies(self.P_cc, self.OF, self.expratio,0,0)
        #h_cc = enthalpies[0]

        #print(h_cc, final_enthalpy)

        #constant specific heats
        #m^2/s^2 = J/kg --> divide by 2000 because 2*1000 to convert
        #print(h_exit, self.C.get_Chamber_Cp(self.P_cc, self.OF, self.expratio,0)*T_exit +((v_exit)**2)/2000, self.C.get_Chamber_Cp(self.P_cc, self.OF, self.expratio)*T_exit +(self.v_exit**2)/2000,( (self.C.get_Chamber_Cp(self.P_cc, self.OF, self.expratio,0)*T_exit +((v_exit)**2)/2000 )- (self.C.get_Chamber_Cp(self.P_cc, self.OF, self.expratio)*T_exit +(self.v_exit**2)/2000) )/ (self.C.get_Chamber_Cp(self.P_cc, self.OF, self.expratio)*T_exit +(self.v_exit**2)/2000))
        #print(self.C.get_Chamber_Cp(self.P_cc, self.OF, self.expratio,0))


        #self.v_exit = (self.v_exit+ v_exit)/2
        
        #print(self.v_exit)


        
        self.instThrust = (self.m_dot_cc_t * self.v_exit) + self.A_exit * (P_exit - self.P_atm)


        #print(self.A_exit * (P_exit - self.P_atm))
        #print((self.m_dot_cc_t * self.v_exit) , self.A_exit * (P_exit - self.P_atm))

        # Recalculate new fuel grain size
        self.m_fuel_t = self.m_fuel_t - self.m_dot_fuel * self.timestep

        self.total_propellant = self.total_propellant + self.m_dot_cc_t * self.timestep

        self.radius = self.radius + self.r_dot_t * self.timestep
        self.A_port_t = np.pi * self.radius**2


        #print(self.m_fuel_t, self.radius)
