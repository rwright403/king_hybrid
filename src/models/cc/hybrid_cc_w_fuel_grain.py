"""

import numpy as np
from src.models.cc._base import BaseChamber
from scipy.optimize import root_scalar
from src.utils.numerical_methods import rk4_step

R_UNIV = 8314  # kJ/mol-K

def partial_derivatives_gamma(C, P_cc, OF, expratio, dP=1e3, dOF=1e-3):
    """
    #Estimate dgamma/dP_cc and dgamma/dOF using finite differences.
    #
    #Parameters:
    #    C      : CEA_Obj or similar equilibrium solver interface
    #    P_cc     : chamber pressure [Pa]
    #    OF       : oxidizer/fuel ratio
    #    expratio : nozzle expansion ratio (needed by CEA)
    #    dP       : perturbation in pressure [Pa]
    #    dOF      : perturbation in mixture ratio
    #Returns:
    #    dgamma_dPcc, dgamma_dOF
"""

    # --- gamma(P_cc, OF) baseline
    #_, gamma0 = cea.get_Chamber_MolWt_gamma(P_cc, OF, expratio)

    # --- dgamma/dP_cc
    _, gamma_plus = C.get_Chamber_MolWt_gamma(P_cc + dP, OF, expratio)
    _, gamma_minus = C.get_Chamber_MolWt_gamma(P_cc - dP, OF, expratio)
    dgamma_dPcc = (gamma_plus - gamma_minus) / (2 * dP)

    # --- dgamma/dOF
    _, gamma_plus = C.get_Chamber_MolWt_gamma(P_cc, OF + dOF, expratio)
    _, gamma_minus = C.get_Chamber_MolWt_gamma(P_cc, OF - dOF, expratio)
    dgamma_dOF = (gamma_plus - gamma_minus) / (2 * dOF)

    return dgamma_dPcc, dgamma_dOF





class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, V_pre_post_cc, m_fuel_i, rho_fuel, a, n, L, A_port,
                 nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        self.C = C

        # fuel + regression
        self.rho_fuel = rho_fuel
        self.a = a
        self.n = n

        # geometry
        self.L = L
        self.r = np.sqrt(A_port / np.pi)
        self.A_port = A_port

        self.A_fuel_grain_od = (m_fuel_i/(rho_fuel*L))+A_port


        self.V_pre_post_cc = V_pre_post_cc

        # initial propellant masses in CC
        self.m_fuel = 0.0
        self.m_ox = 0.0

        # initial values
        self.P_cc = P_atm
        self.P_atm = P_atm



        self.timestep = timestep

        


    def cc_ode_system(self, t, y, m_dot_ox_in):
        
        r, m_ox, m_fuel, P_cc = y #unpack integration var
        print("integration var: ", r, m_ox, m_fuel, P_cc)

        #m_dot_ox_in*=0.05

        ### Solve Thrust:

### BUG!!! PROBLEM IS WITH OF = m_ox/m_fuel AND SOLVER
### WE CAN GET A NEGATIVE OF!!!
        if m_fuel != 0:
            OF = m_ox/m_fuel
        else:
            OF = 2
            m_fuel = 0.0001

        print("OF: ", OF, m_ox, m_fuel)

        

        T_cc = self.C.get_Tcomb(P_cc, OF)
        MW, gamma = self.C.get_Chamber_MolWt_gamma(P_cc, OF, self.nozzle.expratio)

        R_spec = R_UNIV/MW

        print("enter nozzle: ", P_cc, T_cc, gamma, R_spec, OF)
        instThrust, m_dot_exit = self.nozzle.sol_thrust(P_cc, T_cc, gamma, R_spec)

        ### Setup and solve differential eqns:

        # Sol r_dot_ox ODE based on stanford empirical eqn
        A_port = np.pi*r**2
        G_ox = m_dot_ox_in/A_port
        r_dot =self.a*G_ox**self.n


        #print("OF", OF, m_fuel, m_ox, G_ox)


        # Solve m_dot_fuel, m_dot_ox

        V_dot = (2*np.pi*r*self.L)*r_dot
        m_dot_fuel_in = self.rho_fuel * V_dot # Know m_dot_ox_in from tank

        m_dot_fuel = m_dot_fuel_in - (m_dot_exit)/(1+(1/OF))
        m_dot_ox = m_dot_ox_in - (m_dot_exit)/(1+OF)

        print("m_dot_ox: ", m_dot_ox, m_dot_fuel)

        # Solve P_dot_cc
        #m_dot_cc = m_dot_fuel + m_dot_ox
        m_dot_cc = m_dot_ox_in + m_dot_fuel_in - m_dot_exit
        print("m_dot_cc! ", m_dot_cc, (m_dot_ox_in+m_dot_fuel_in) , m_dot_exit)

        V_cc = self.L*(self.A_fuel_grain_od - A_port) + self.V_pre_post_cc
        #print("V_cc sign: ", V_cc)


        dgamma_dPcc, dgamma_dOF = partial_derivatives_gamma(self.C, P_cc, OF, self.nozzle.expratio)
        
        #NOTE: CHECK UNITS ON CP, I THINK ITS IN KJ/KG K RN
        cp = self.C.get_Chamber_Cp(P_cc, OF, self.nozzle.expratio)*1000
        print("cp inputs: ", cp, P_cc, OF, self.nozzle.expratio)

        #print("entering P_dot eqn, check denoms: ",(1 - (P_cc/(V_cc-1))*dgamma_dPcc), V_cc, ((V_cc-1)*m_fuel))

        P_dot =  1/(1 - (P_cc/(V_cc-1)) * dgamma_dPcc) * ( ((gamma-1)/V_cc)*m_dot_cc*cp*T_cc - gamma*(P_cc/V_cc)*V_dot + (P_cc/((V_cc-1)*m_fuel)) * dgamma_dOF*(m_dot_ox - OF*m_dot_fuel) )

        #print("P_cc: ", P_cc)
        print("P_dot: ", P_dot,  ((gamma-1)/V_cc)*m_dot_cc*cp*T_cc, - gamma*(P_cc/V_cc)*V_dot , + (P_cc/((V_cc-1)*m_fuel))*dgamma_dOF*(m_dot_ox-OF*m_dot_fuel) )
        #print("P_dot first term: ",((gamma-1)/V_cc)*m_dot_cc*cp*T_cc, ((gamma-1)/V_cc), m_dot_cc,cp,T_cc)
        #print("P_dot second_term: ", -gamma*(P_cc/V_cc)*V_dot, gamma, P_cc,V_cc*V_dot)

        #print("m_dot: ", m_dot_cc, m_dot_exit, (m_ox+m_fuel))



        return [r_dot, m_dot_ox, m_dot_fuel, P_dot], {"P_cc": P_cc, "thrust": instThrust}
    

    def cc_ode_system_rk(self, t, y, m_dot_ox):
        out, _ = self.cc_ode_system(t, y, m_dot_ox)
        return out



    def inst(self, m_dot_ox, _: float = 0.0): #NOTE:m_dot_fuel solved from fuel grain in this cc. This is an artifact so program is modular
        
        t = 0
        y0 = [self.r, self.m_ox, self.m_fuel, self.P_cc]
        y_new = rk4_step(self.cc_ode_system_rk, 0.0, y0, self.timestep, m_dot_ox)
        self.r, self.m_ox, self.m_fuel, self.P_cc = y_new

        #print("fuel/ox mass: ", self.m_fuel, self.m_ox)

        _, out = self.cc_ode_system(t, y_new, m_dot_ox)
        return out    
    

"""



import numpy as np
from src.models.cc._base import BaseChamber
from src.utils.numerical_methods import rk4_step

R_UNIV = 8314  # J/kmol-K

def partial_derivatives_gamma(C, P_cc, OF, expratio, dP=1e3, dOF=1e-3):
    #Estimate dgamma/dP_cc and dgamma/dOF using finite differences.
    _, gamma_plus = C.get_Chamber_MolWt_gamma(P_cc + dP, OF, expratio)
    _, gamma_minus = C.get_Chamber_MolWt_gamma(P_cc - dP, OF, expratio)
    dgamma_dPcc = (gamma_plus - gamma_minus) / (2 * dP)

    _, gamma_plus = C.get_Chamber_MolWt_gamma(P_cc, OF + dOF, expratio)
    _, gamma_minus = C.get_Chamber_MolWt_gamma(P_cc, OF - dOF, expratio)
    dgamma_dOF = (gamma_plus - gamma_minus) / (2 * dOF)

    return dgamma_dPcc, dgamma_dOF


class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, V_pre_post_cc, m_fuel_i, rho_fuel, a, n, L, A_port,
                 nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        self.C = C

        # fuel regression parameters
        self.rho_fuel = rho_fuel
        self.a = a
        self.n = n

        # geometry
        self.L = L
        self.r = np.sqrt(A_port / np.pi)
        self.A_port = A_port
        self.A_fuel_grain_od = (m_fuel_i / (rho_fuel * L)) + A_port
        self.r_fuel_grain_outer = np.sqrt(self.A_fuel_grain_od/np.pi)
        self.V_pre_post_cc = V_pre_post_cc

        # chamber state variables
        self.m_cc = 0.0       # total mass in chamber (fuel + ox)
        self.V_cc = 0.0       # initialized here, will be sol for
        self.P_cc = P_atm     # start at ambient
        self.P_atm = P_atm
        self.timestep = timestep

#NOTE: THIS IS A QUICK MODEL, I DID NOT TRACK FUEL GRAIN MASS, SO CASES WHERE FUEL GRAIN FULLY BURNS, THIS MODEL WILL FAIL

    def cc_ode_system(self, t, y, m_dot_ox_in):
        r, m_cc, P_cc = y

        # --- Fuel inflow from grain regression
        A_port = np.pi * r**2
        G_ox = m_dot_ox_in / A_port
        r_dot = self.a * G_ox**self.n

        V_dot = (2 * np.pi * r * self.L) * r_dot
        m_dot_fuel_in = self.rho_fuel * V_dot

        # --- Mixture ratio from inflows only
        OF = m_dot_ox_in / max(m_dot_fuel_in, 1e-6)
        #print("OF", OF, G_ox)

        # --- Thermochem
        T_cc = self.C.get_Tcomb(P_cc, OF)
        MW, gamma = self.C.get_Chamber_MolWt_gamma(P_cc, OF, self.nozzle.expratio)
        R_spec = R_UNIV / MW

        instThrust, m_dot_exit = self.nozzle.sol_thrust(P_cc, T_cc, gamma, R_spec)
        

        # --- Chamber mass balance (total)
        m_dot_cc = m_dot_ox_in + m_dot_fuel_in - m_dot_exit

        print("m_dot_cc: ", m_dot_cc, m_dot_ox_in, m_dot_fuel_in, - m_dot_exit )

        # --- Chamber volume
        self.V_cc = self.L * A_port + self.V_pre_post_cc
        


        print("V_cc: ", self.V_cc, (self.A_fuel_grain_od - A_port), self.A_fuel_grain_od, - A_port)


        # --- Gamma derivatives
        dgamma_dPcc, dgamma_dOF = partial_derivatives_gamma(self.C, P_cc, OF, self.nozzle.expratio)

        #print("partials: ", dgamma_dPcc, dgamma_dOF)

        # --- cp from CEA (convert kJ/kg-K -> J/kg-K)
        cp = self.C.get_Chamber_Cp(P_cc, OF, self.nozzle.expratio) * 1000

        #m_fuel_chamber = m_cc/OF
        #NOTE: added this before it was just m_cc, also this is just in cc, not the fuel grain accounting!

        # --- Pressure ODE (Zimmerman/McGill style form)
        P_dot = 1/(1 - (P_cc /(gamma-1))*dgamma_dPcc) * ( ((gamma-1)/self.V_cc)*m_dot_cc*cp*T_cc - gamma*(P_cc/self.V_cc)*V_dot )#+ (P_cc/((gamma-1)*max(m_fuel, 1e-6)))*dgamma_dOF*(m_dot_ox_in-OF*m_dot_fuel_in) )


        #print("vars: ", m_dot_exit, OF, P_dot, ((gamma-1)/V_cc)*m_dot_cc*cp*T_cc, - gamma*(P_cc/V_cc)*V_dot , (P_cc/((V_cc-1)*max(m_fuel_chamber, 1e-6)))*dgamma_dOF*(m_dot_ox_in-OF*m_dot_fuel_in) )
        
        #NOTE: INSTABILITY IN THIS TERM: ((gamma-1)/V_cc)*m_dot_cc*cp*T_cc
        #print("m_fuel_chamber: ", m_fuel_chamber, m_dot_cc, cp, T_cc)

        #(P_cc / ((V_cc - 1) * m_dot_fuel_in)) * dgamma_dOF * (m_dot_ox_in - OF*m_dot_fuel_in)

        #print("P_dot: ", P_dot,  ((gamma-1)/V_cc)*m_dot_cc*cp*T_cc, - gamma*(P_cc/V_cc)*V_dot , + (P_cc/((V_cc-1)*max(m_fuel, 1e-6)))*dgamma_dOF*(m_dot_ox_in-OF*m_dot_fuel_in) )
        


        return [r_dot, m_dot_cc, P_dot], {"P_cc": P_cc, "thrust": instThrust}

    def cc_ode_system_rk(self, t, y, m_dot_ox):
        out, _ = self.cc_ode_system(t, y, m_dot_ox)
        return out

    def inst(self, m_dot_ox, _: float = 0.0):
        y0 = [self.r, self.m_cc, self.P_cc]
        y_new = rk4_step(self.cc_ode_system_rk, 0.0, y0, self.timestep, m_dot_ox)
        self.r, self.m_cc, self.P_cc = y_new

        _, out = self.cc_ode_system(0, y_new, m_dot_ox)

        print("P_cc: ", self.P_cc, self.m_cc)
        return out
    


"""

import numpy as np
from src.models.cc._base import BaseChamber
from src.utils.numerical_methods import rk4_step

R_UNIV = 8314  # J/kmol-K

class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, V_pre_post_cc, m_fuel_i, rho_fuel, a, n, L, A_port,
                 nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        self.C = C
        self.rho_fuel = rho_fuel
        self.a = a
        self.n = n
        self.L = L
        self.r = np.sqrt(A_port / np.pi)
        self.A_port = A_port
        self.A_fuel_grain_od = (m_fuel_i / (rho_fuel * L)) + A_port
        self.V_pre_post_cc = V_pre_post_cc

        # chamber states
        self.m_ox_c = 0.0
        self.m_fuel_c = 0.0
        self.P_cc = P_atm
        self.gamma_prev = 1.2  # safe initial guess
        self.timestep = timestep

    def cc_ode_system(self, t, y, m_dot_ox_in):
        r, m_ox_c, m_fuel_c, P_cc = y

        # --- regression rate
        A_port = np.pi * r**2
        G_ox = m_dot_ox_in / A_port
        r_dot = self.a * G_ox**self.n

        V_dot = (2*np.pi*r*self.L) * r_dot
        m_dot_fuel_in = self.rho_fuel * V_dot

        # OF ratio
        OF = m_dot_ox_in / max(m_dot_fuel_in, 1e-6)

        # thermochem
        T_cc = self.C.get_Tcomb(P_cc, OF)
        MW, gamma = self.C.get_Chamber_MolWt_gamma(P_cc, OF, self.nozzle.expratio)
        R_spec = R_UNIV / MW

        instThrust, m_dot_exit = self.nozzle.sol_thrust(P_cc, T_cc, gamma, R_spec)

        # --- chamber mass balances
        m_dot_ox_out = (m_dot_exit) / (1 + OF)
        m_dot_fuel_out = (m_dot_exit) / (1 + 1/OF)

        dMox_dt = m_dot_ox_in - m_dot_ox_out
        dMfuel_dt = m_dot_fuel_in - m_dot_fuel_out

        # --- chamber volume
        V_cc = self.L * (self.A_fuel_grain_od - A_port) + self.V_pre_post_cc

        # --- gamma time derivative (Bernarde uses one-sided differencing)
        dgamma_dt = (gamma - self.gamma_prev) / self.timestep
        self.gamma_prev = gamma

        cp = self.C.get_Chamber_Cp(P_cc, OF, self.nozzle.expratio) * 1000 # convert units from [kJ/(kg K)] to [J/(kg K)]

        # --- chamber pressure ODE (Bernarde eqn 17)
        m_dot_cc = m_dot_ox_in + m_dot_fuel_in - m_dot_exit
        P_dot = ((gamma - 1)/V_cc) * (m_dot_cc * cp * T_cc) \
                - (gamma * P_cc / V_cc) * V_dot \
                + (P_cc / (V_cc - 1)) * dgamma_dt

        return [r_dot, dMox_dt, dMfuel_dt, P_dot], {"P_cc": P_cc, "thrust": instThrust}

    def cc_ode_system_rk(self, t, y, m_dot_ox):
        out, _ = self.cc_ode_system(t, y, m_dot_ox)
        return out

    def inst(self, m_dot_ox, _: float = 0.0):
        y0 = [self.r, self.m_ox_c, self.m_fuel_c, self.P_cc]
        y_new = rk4_step(self.cc_ode_system_rk, 0.0, y0, self.timestep, m_dot_ox)
        self.r, self.m_ox_c, self.m_fuel_c, self.P_cc = y_new

        print("P_cc: ", self.P_cc, (self.m_ox_c+self.m_fuel_c) )

        _, out = self.cc_ode_system(0, y_new, m_dot_ox)
        return out
"""