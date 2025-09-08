import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
from src.models.cc._base import BaseChamber


#TODO: DOUBLE CHECK PropsSI units  



class adiabatic_lre_cc_model(BaseChamber):

    def __init__(self, oxidizer_name, fuel_name, nozzle, P_atm, C, timestep):
        self.nozzle = nozzle

        #use SI UNITS
        self.C = CEA_Obj(oxName=oxidizer_name, fuelName=fuel_name, pressure_units='Pa', isp_units='sec', cstar_units='m/s', temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='kJ/kg', density_units='kg/m^3', specific_heat_units='kJ/kg-K')
        #NOTE: SEEMS TO BREAK ABOVE

        self.P_cc = 0
        self.OF = 0

        self.y = 0 
        self.R = 0 # J/kg*K

        # initial values (treat as private members)
        self.P_cc = P_atm
        self.P_atm = P_atm
        self.v_exit = 0
        self.r_dot_t = 0.1 

        self.m_dot_cc_t = 0
        self.prev_thrust = 0
        self.instThrust = 0

        self.A_throat = nozzle.A_throat 
        self.A_exit = nozzle.A_exit 
        self.expratio = self.A_exit / self.A_throat

        self.TIMESTEP = timestep


    def inst(self, m_dot_ox, m_dot_fuel):
        
        #print("this should not be nan: ",m_dot_ox,m_dot_fuel)
        #update O/F ratio and m_dot_cc
        self.OF = m_dot_ox/m_dot_fuel
        self.m_dot_cc_t = m_dot_ox + m_dot_fuel
        #print("testing with a factor multiplied to m_dot_cc_t just to feel")
        #self.m_dot_cc_t *= 2
        
        #CALL CEA to solve combustion
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K)
        self.y = fluid_prop[1] # (-)
        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.expratio, 0, 1)
        T_cc = temperatures[0]

        ##("TEMPERATURE: ",T_cc)
        
        ##(self.m_dot_cc_t, m_dot_fuel, m_dot_ox)

        #NOTE: P_cc too high? --> too low
        #print("AAAAA:", self.m_dot_cc_t, self.A_throat, self.R, T_cc, self.y)
        self.P_cc = (self.m_dot_cc_t / self.A_throat ) * np.sqrt( self.R*T_cc )  / ( np.sqrt( self.y * (2 / (self.y+1))**( (self.y+1)/(self.y-1) ) ) )

        #print("NOTE: for debugging multiplied P_cc by 2")

        #print("in cc:", self.P_cc)

        #resolve fluid properties
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K)
        self.y = fluid_prop[1] # (-)

        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.expratio, 0, 1)
        T_cc = temperatures[0]
        #T_throat = temperatures[1]
        #T_exit = temperatures[2]
        #print(T_cc, T_throat, T_exit)


        ###START: NOZZLE
        thrust = self.nozzle.compute(
                    P_cc=self.P_cc,
                    T_cc=T_cc,
                    m_dot=self.m_dot_cc_t,
                    gamma=self.y,
                    R=self.R,
                    chamber_obj=self.C
        )
        
        #solve thrust
        self.prev_thrust = self.instThrust
        self.instThrust = thrust

        return {"P_cc": self.P_cc, "thrust": self.instThrust}