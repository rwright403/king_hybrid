#https://rocketcea.readthedocs.io/en/latest/functions.html#module-rocketcea.cea_obj
import numpy as np
from scipy.optimize import root_scalar
from rocketcea.cea_obj import add_new_fuel
from rocketcea.cea_obj_w_units import CEA_Obj
from src.models.cc._base import BaseChamber
from src.utils.numerical_methods import rk4_step

R_UNIV = 8314  # J / (kmol K)




class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, L_star, m_fuel_i, rho_fuel, a, n, L, A_port, nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        self.C = C

        # fuel grain properties
        self.rho_fuel = rho_fuel
        self.a = a
        self.n = n
        #self.h_latent_fuel = 400e3  # J/kg, placeholder

        # fuel grain geometry and mass
        self.m_fuel_t = m_fuel_i
        self.L = L
        #self.A_port = A_port
        self.r = np.sqrt(A_port / np.pi)

        print("r!", self.r)


        # oxidizer CEA object
        ox_card = """#fuel paraffin  N 2   O 1    wt%=100.00"""
        add_new_fuel('N2O', ox_card)
        self.oxCEA = CEA_Obj(oxName='N2O', fuelName='N2O')

        self.timestep = timestep
        self.m_dot_fuel = 0.0





    def inst(self, m_dot_ox, _):
        # --- 1. Regression
        A_port = np.pi * self.r**2
        G_ox = m_dot_ox / max(A_port, 1e-12)
        print("G_ox: ", G_ox)
        r_dot = self.a * G_ox**self.n
        m_dot_fuel = self.rho_fuel * (2*np.pi*self.r*self.L) * r_dot
        OF = m_dot_ox / max(m_dot_fuel, 1e-9)

        # --- 2. Residual closure
        def mass_balance(P_cc):
            T_cc = self.C.get_Tcomb(P_cc, OF)
            MW, gamma = self.C.get_Chamber_MolWt_gamma(P_cc, OF,
                                                       self.nozzle.expratio)
            R_spec = R_UNIV / MW

            #print("thermochemistry: ", P_cc, T_c, gamma, R_spec)
            _, m_dot_noz = self.nozzle.sol_thrust(P_cc, T_cc, gamma, R_spec)

            #print("res: ", ((m_dot_ox + m_dot_fuel) - m_dot_noz), m_dot_ox, m_dot_fuel, P_cc, OF)
            return (m_dot_ox + m_dot_fuel) - m_dot_noz

        sol = root_scalar(mass_balance,
                          bracket=[1e5, 1e8],
                          method="bisect")
        if not sol.converged:
            raise RuntimeError("Pcc solve failed")
        self.P_cc = sol.root

        # --- 3. Final thermo + thrust
        T_cc = self.C.get_Tcomb(self.P_cc, OF)
        MW, gamma = self.C.get_Chamber_MolWt_gamma(self.P_cc, OF,
                                                   self.nozzle.expratio)
        R_spec = R_UNIV / MW
        self.instThrust, _ = self.nozzle.sol_thrust(self.P_cc, T_cc, gamma, R_spec)

        # --- 4. Update port radius
        self.r += r_dot * self.timestep

        return {"P_cc": self.P_cc, "thrust": self.instThrust}
    
