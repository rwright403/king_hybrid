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

        #print("m_dot_cc: ", m_dot_cc, m_dot_ox_in, m_dot_fuel_in, - m_dot_exit )

        # --- Chamber volume
        self.V_cc = self.L * A_port + self.V_pre_post_cc
        


        #print("V_cc: ", self.V_cc, (self.A_fuel_grain_od - A_port), self.A_fuel_grain_od, - A_port)


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
        


        return [r_dot, m_dot_cc, P_dot], {"P_cc": P_cc, "thrust": instThrust, "m_dot_fuel": m_dot_fuel_in}

    def cc_ode_system_rk(self, t, y, m_dot_ox):
        out, _ = self.cc_ode_system(t, y, m_dot_ox)
        return out

    def inst(self, m_dot_ox, _: float = 0.0):
        y0 = [self.r, self.m_cc, self.P_cc]
        y_new = rk4_step(self.cc_ode_system_rk, 0.0, y0, self.timestep, m_dot_ox)
        self.r, self.m_cc, self.P_cc = y_new

        print(f"             |  cc: {self.r:.3f}, {self.m_cc:.3f}, {self.P_cc:.0f}")

        _, out = self.cc_ode_system(0, y_new, m_dot_ox)

        return out
