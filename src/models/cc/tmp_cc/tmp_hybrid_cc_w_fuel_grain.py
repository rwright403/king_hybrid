#https://rocketcea.readthedocs.io/en/latest/functions.html#module-rocketcea.cea_obj
import numpy as np
from scipy.optimize import root_scalar
from rocketcea.cea_obj import add_new_fuel
from rocketcea.cea_obj_w_units import CEA_Obj
from src.models.cc._base import BaseChamber
from src.utils.numerical_methods import rk4_step

R_UNIV = 8314 # J / (kg K)

def cea_props(P_cc, C, OF, expratio):
    T = C.get_Tcomb(P_cc, OF)
    MW, gamma = C.get_Chamber_MolWt_gamma(P_cc, OF, expratio)
    R_spec = (R_UNIV/MW)
    h = C.get_Chamber_H(P_cc, OF, expratio)
    return T, R_spec, gamma, h


def uerror_cea(P_cc, C, OF, expratio, u):
    T, R_spec, _, h = cea_props(P_cc, C, OF, expratio)
    print("cea props: ", T, R_spec, h )
    u_cea = h - R_spec*T
    print("sol u_cea: ", h, - R_spec*T)
    print("uerr: ", u_cea - u, u_cea, u)
    return u_cea - u

class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, L_star, m_fuel_i, rho_fuel, a, n, L, A_port, nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        #use SI UNITS
        self.C = C

        # fuel grain material
        self.rho_fuel = rho_fuel 
        self.a = a 
        self.n = n 
        self.h_latent_fuel = 400e3 #NOTE: TMP

        #fuel grain geometry and mass
        self.m_fuel_t = m_fuel_i 
        self.L = L
        self.A_port = A_port
        self.r_prev = np.sqrt(A_port / np.pi)


        # initial values
        self.P_cc = P_atm
        self.P_atm = P_atm
        self.v_exit = 0
        self.r_dot_t = 0.1 

        self.m_dot_reacs = 0   

        self.V_cc = L_star * self.nozzle.A_throat

        print("init check vol: ", self.V_cc, L_star, self.nozzle.A_throat)

        self.timestep = timestep    

        # setup
        self.m_dot_fuel = 0

        self.m_cc = 1e-6
        T, R_spec, _, h = cea_props(self.P_cc, C, 0.75, self.nozzle.expratio)
        self.U_cc = self.m_cc*(h - R_spec*T)

        ## setup oxidizer CEA obj for enthalpy consistency:
        ox_card = """#fuel paraffin  N 2   O 1    wt%=100.00"""
        add_new_fuel('N2O', ox_card)
        self.oxCEA = CEA_Obj(oxName='N2O', fuelName='N2O')

    def cc_ode_system(self, t, y, m_dot_ox):
        
        #TODO: unpack rk var
        m_cc, U_cc, r = y

        #sol r_dot_ox
        A_port = np.pi*r**2
        G_ox = m_dot_ox/A_port
        r_dot =self.a*G_ox**self.n

        #solve m_dot_fuel to solve O/F
        m_dot_fuel = self.rho_fuel*(2*np.pi*r*self.L)*r_dot
        OF = m_dot_ox / max(m_dot_fuel, 1e-9) #on startup m_dot_ox = 0

        #solve T_cc
        u_cc = U_cc/max(m_cc, 1e-3) #on startup m_cc = 0)
        print("u_cc", u_cc)

        sol = root_scalar(
            lambda P_guess: uerror_cea(P_guess, self.C, OF, self.nozzle.expratio, u_cc),
            method="secant",
            x0=self.P_cc,
            x1=self.P_cc*1.1,
            xtol=1000,
            maxiter=1000
        )
        if not sol.converged:
            raise RuntimeError("uerror_cea solver did not converge")
        self.P_cc = sol.root

        print("P_cc from uerror_cea: ", self.P_cc)

        # call CEA for thermo properties:
        T_cc, R_spec, gamma, h_cc = cea_props(self.P_cc, self.C, OF, self.nozzle.expratio)

        # sol m_dot out of cc from nozzler
        instThrust, m_dot_reacs = self.nozzle.sol_thrust(self.P_cc, T_cc, gamma, R_spec)

        # mass and energy balance:
        m_dot_cc = m_dot_fuel + m_dot_ox - m_dot_reacs
        
#BUG: ENTHALPY REFERENCE MUST BE CONSISTENT BETWEEN CEA AND NITROUS TERMS HERE
        #U_dot_cc = m_dot_ox*(h_ox + (v_ox**2)/2) +m_dot_fuel*h_fuel - m_dot_reacs*(h_cea + (v_reactants_out**2)/2)
        #NOTE: neglecting velocity term of RTT for simplicity
        h_ox = self.oxCEA.get_Chamber_H(self.P_cc, OF, self.nozzle.expratio)

        U_dot_cc = m_dot_ox*h_ox +m_dot_fuel*self.h_latent_fuel - m_dot_reacs*h_cc

        return [m_dot_cc, U_dot_cc, r_dot], {"P_cc": self.P_cc, "thrust": instThrust}
    

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







from src.models._thermo.n2o_thermo_span_wagner_class import SpanWagnerEOS_SingleState
def compute_reference_offset(C,sweos, T_ref, OF=6.0, expratio=40.0):
    """Compute delta_u = u_SW - u_CEA at reference state."""

    P_ref = sweos.P

    # --- CEA values
    MW, gamma = C.get_Chamber_MolWt_gamma(P_ref, OF, expratio)
    R_spec = R_UNIV / MW
    h_cea = C.get_Chamber_H(P_ref, OF, expratio)
    u_cea = h_cea - R_spec * T_ref

    print("difference in h: ", (h_cea-sweos.h), h_cea, sweos.h)

    # --- Span–Wagner values (NIST convention)
    u_sw = sweos.u  # J/kg

    delta_u = u_sw - u_cea
    return delta_u


ox_card = """#fuel paraffin  N 2   O 1    wt%=100.00"""
add_new_fuel('N2O', ox_card)
oxCEA = CEA_Obj(oxName='N2O', fuelName='N2O')

T_ref = 1000
sweos = SpanWagnerEOS_SingleState(100, T_ref)


print("ref int energy: ", compute_reference_offset(oxCEA,sweos, T_ref) )



"""
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
#BUG: THIS SEEMS WRONG BELOWok       
            self.m_dot_fuel = self.rho_fuel * np.pi * ((self.r_dot_t + np.sqrt(self.A_port_t / np.pi))**2 - (self.A_port_t / np.pi)) * self.L
            
            if self.m_fuel_t <= 0:
                self.m_dot_fuel = 0
                self.OF = 0
            else:
                self.OF = m_dot_ox / self.m_dot_fuel

            i += 1

            #fixed point iteration on P_cc

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

            print("nozzle in: ", self.P_cc, T_cc, y, R_spec )

            instThrust, self.m_dot_reactants_out = self.nozzle.sol_thrust(self.P_cc, T_cc, y, R_spec)

            # As it loops it will update mass balance

        ### Recalculate new fuel grain size for next timestep
        self.m_fuel_t = self.m_fuel_t - self.m_dot_fuel * self.timestep

        self.radius = self.radius + self.r_dot_t * self.timestep
        self.A_port_t = np.pi * self.radius**2

        return {"P_cc": self.P_cc, "thrust": instThrust}
"""
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
    def __init__(self, oxidizer_name, fuel_name, fuel_properties, m_fuel_i, rho_fuel, a, n, L, A_port, P_atm, A_throat, A_exit, timestep):
       
"""