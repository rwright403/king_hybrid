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


class adiabatic_lre_cc_model(BaseChamber):

    def __init__(self, V_cc, nozzle, P_atm, C, timestep):

        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        self.C = C

        # geometry
        self.V_cc = V_cc #NOTE: use L_star or V_cc?

        # chamber state variables
        self.m_cc = 0.0       # total mass in chamber (fuel + ox)
        self.P_cc = P_atm     # start at ambient
        self.P_atm = P_atm #NOTE: unused? delete?
        self.timestep = timestep

    def cc_ode_system(self, t, y, m_dot_ox_in, m_dot_fuel_in):
        m_cc, P_cc = y

        # --- Mixture ratio from inflows only
        OF = m_dot_ox_in / max(m_dot_fuel_in, 1e-6)
        #print("OF", OF, G_ox)

        # --- Thermochem

        #print("enter thermochem: ", P_cc, OF, m_dot_ox_in, m_dot_fuel_in, m_cc)

        T_cc = self.C.get_Tcomb(P_cc, OF)
        MW, gamma = self.C.get_Chamber_MolWt_gamma(P_cc, OF, self.nozzle.expratio)
        R_spec = R_UNIV / MW

        instThrust, m_dot_exit = self.nozzle.sol_thrust(P_cc, T_cc, gamma, R_spec)
        

        # --- Chamber mass balance (total)
        m_dot_cc = m_dot_ox_in + m_dot_fuel_in - m_dot_exit


        # --- Gamma derivatives
        dgamma_dPcc, _ = partial_derivatives_gamma(self.C, P_cc, OF, self.nozzle.expratio)

        #print("partials: ", dgamma_dPcc, dgamma_dOF)

        # --- cp from CEA (convert kJ/kg-K -> J/kg-K)
        cp = self.C.get_Chamber_Cp(P_cc, OF, self.nozzle.expratio) * 1000

        #m_fuel = m_cc/OF
        #NOTE: added this before it was just m_cc

        # --- Pressure ODE (Zimmerman/McGill style form)
        P_dot = 1/(1 - (P_cc /(gamma-1))*dgamma_dPcc) * ( ((gamma-1)/self.V_cc)*m_dot_cc*cp*T_cc ) #+ (P_cc/((V_cc-1)*max(m_fuel, 1e-6)))*dgamma_dOF*(m_dot_ox_in-OF*m_dot_fuel_in) )


        #print("vars: ", m_dot_exit, OF, P_dot, ((gamma-1)/self.V_cc)*m_dot_cc*cp*T_cc, (P_cc/((self.V_cc-1)*max(m_fuel, 1e-6)))*dgamma_dOF*(m_dot_ox_in-OF*m_dot_fuel_in) )


        return [m_dot_cc, P_dot], {"P_cc": P_cc, "thrust": instThrust}

    def cc_ode_system_rk(self, t, y, m_dot_ox, m_dot_fuel):
        out, _ = self.cc_ode_system(t, y, m_dot_ox, m_dot_fuel)
        return out

    def inst(self, m_dot_ox, m_dot_fuel):
        y0 = [self.m_cc, self.P_cc]
        y_new = rk4_step(self.cc_ode_system_rk, 0.0, y0, self.timestep, m_dot_ox, m_dot_fuel)
        self.m_cc, self.P_cc = y_new

        print(f"             |  cc: {self.m_cc:.3f}, {self.P_cc:.3f}")

        _, out = self.cc_ode_system(0, y_new, m_dot_ox, m_dot_fuel)

        print("P_cc: ", self.P_cc, self.m_cc)
        return out


"""
class adiabatic_lre_cc_model(BaseChamber):

    def __init__(self, L_star, nozzle, P_atm, C, timestep):
        self.nozzle = nozzle

        #use SI UNITS
        self.C = C

        self.OF = 0.5

        # initial values
        self.P_atm = P_atm
        self.P_cc = P_atm
        self.v_exit = 0
        self.r_dot_t = 0.1 

        self.m_cc = 0
        self.m_dot_reactants_out = 0

        self.V_cc = L_star * self.nozzle.A_throat

        self.timestep = timestep


    def inst(self, m_dot_ox, m_dot_fuel):

        #fixed point iteration on P_cc

        for _ in range(10):

            #Update mass balance: 
            OF = m_dot_ox / m_dot_fuel
            m_dot_propellants_in = m_dot_ox + m_dot_fuel

            self.m_cc += (m_dot_propellants_in-self.m_dot_reactants_out)*self.timestep

            # CEA to solve combustion properties
            MW, y = self.C.get_Chamber_MolWt_gamma(self.P_cc, OF, self.nozzle.expratio)
            T_cc = self.C.get_Tcomb(self.P_cc, OF)

            # Treat reactant fluid as ideal gas and resolve P_cc
            R_spec = (R_UNIV/MW)
            self.P_cc = R_spec * (T_cc*self.m_cc/self.V_cc)

            instThrust, self.m_dot_reactants_out = self.nozzle.sol_thrust(self.P_cc, T_cc, y, R_spec)

            # As it loops it will update mass balance and converge
            #print(self.P_cc, instThrust)

        return {"P_cc": self.P_cc, "thrust": instThrust}
"""