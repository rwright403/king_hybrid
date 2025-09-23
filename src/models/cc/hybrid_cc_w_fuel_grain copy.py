import numpy as np
from src.models.cc._base import BaseChamber
from src.models._thermo.cantera_lookup import *
from scipy.optimize import root_scalar
from src.utils.numerical_methods import rk4_step

R_UNIV = 8.314  # J/mol-K

class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, L_star, m_fuel_i, rho_fuel, a, n, L, A_port,
                 nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        self.thermochemistry = CombustionChamberLookup(C) #TODO: USING PLACEHOLDERS IN cantera_lookup.py , need to refactor and replace 




        # fuel + regression
        self.rho_fuel = rho_fuel
        self.a = a
        self.n = n

        # geometry
        self.L = L
        self.radius = np.sqrt(A_port / np.pi)
        self.A_port_t = A_port
        self.V_cc = L_star * self.nozzle.A_throat

        # masses
        self.m_fuel_t = m_fuel_i
        self.m_cc = 1e-6  # tiny fill to avoid div by zero

        # initial values
        self.P_cc = P_atm
        self.P_atm = P_atm
        #self.m_dot_fuel = 0.0 
        #self.m_dot_reactants_out = 0.0
        #self.OF = 0.0
        self.timestep = timestep

        #self.T_cc = 350


    def cc_ode_system(self, t, y, m_dot_ox):
        
        m_cc, U_cc, r = y #unpack integration var

        #sol r_dot_ox ODE based on stanford empirical eqn
        A_port = np.pi*r**2
        G_ox = m_dot_ox/A_port
        r_dot =self.a*G_ox**self.n

        #solve m_dot_fuel to solve O/F
        m_dot_fuel = self.rho_fuel*(2*np.pi*r*self.L)*r_dot
        OF = m_dot_ox / max(m_dot_fuel, 1e-9) #on startup m_dot_ox = 0


        u = U_cc/m_cc
        rho = m_cc/self.V_cc
        T_cc, P_cc, h_cc, cp, gamma, MW = self.thermochemistry.query(u,rho,OF)

        R_spec = R_UNIV/MW
        instThrust, m_dot_exit = self.nozzle.sol_thrust(P_cc, T_cc, gamma, R_spec)

        #sol mass and energy balance
        m_dot_cc = (m_dot_ox+m_dot_fuel) - m_dot_exit





#NOTE: FOR FIRST PASS NEGLECTING THE VELOCITY TERMS!!!
        U_dot_cc = m_dot_ox*h_ox +m_dot_fuel*self.h_latent_fuel - m_dot_exit*h_cc


        return [m_dot_cc, U_dot_cc, r_dot], {"P_cc": P_cc, "thrust": instThrust}
    

    def cc_ode_system_rk(self, t, y, m_dot_ox):
        out, _ = self.cc_ode_system(t, y, m_dot_ox)
        return out



    def inst(self, m_dot_ox, _: float = 0.0): #NOTE:m_dot_fuel solved from fuel grain in this cc. This is an artifact so program is modular
        
        t = 0
#TODO: MAKE SURE THESE END UP IN INIT
        y0 = [self.m_cc, self.U_cc, self.r_prev]
        y_new = rk4_step(self.cc_ode_system_rk, 0.0, y0, self.timestep, m_dot_ox)
        self.m_cc, self.U_cc, self.r_prev = y_new

        self.A_port = np.pi*self.r_prev**2#NOTE: DO WE NEED THIS?

        _, out = self.cc_ode_system(t, y_new, m_dot_ox)
        return out    