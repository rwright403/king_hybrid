#https://rocketcea.readthedocs.io/en/latest/functions.html#module-rocketcea.cea_obj
import numpy as np
from scipy.optimize import root_scalar
from rocketcea.cea_obj import add_new_fuel
from rocketcea.cea_obj_w_units import CEA_Obj
from src.models.cc._base import BaseChamber
from src.utils.numerical_methods import rk4_step

R_UNIV = 8314  # J / (kmol K)

class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, L_star, m_fuel_i, rho_fuel, a, n, L, A_port,
                 nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        self.C = C

        # fuel grain properties
        self.rho_fuel = rho_fuel
        self.a = a
        self.n = n

        # fuel grain geometry and mass
        self.m_fuel_t = m_fuel_i
        self.L = L
        self.r = np.sqrt(A_port/np.pi)

        # chamber volume
        self.V_cc = L_star * self.nozzle.A_throat

        self.P_cc = P_atm
        self.cv = 1000.0  # J/kg-K, seed guess

        # oxidizer CEA object
        ox_card = """#fuel paraffin  N 2   O 1    wt%=100.00"""
        add_new_fuel('N2O', ox_card)
        self.oxCEA = CEA_Obj(oxName='N2O', fuelName='N2O')

        self.timestep = timestep
        self.m_dot_fuel = 0.0

    
    def inst(self, m_dot_ox, _: float = 0.0): #NOTE:m_dot_fuel solved from fuel grain in this cc. This is an artifact so program is modular

        # Geometry update
        A_port = np.pi *self.r**2
        G_ox = m_dot_ox / max(A_port, 1e-12)
        r_dot = self.a * G_ox**self.n

        # Fuel mass flow
        m_dot_fuel = self.rho_fuel * (2 * np.pi * self.r * self.L) * r_dot
        OF = m_dot_ox / max(m_dot_fuel, 1e-9)

        # function for mass balance residual
        def mass_balance(P_guess):
            # Get c* from CEA
            cstar = self.C.get_Cstar(P_guess, OF)
            mdot_noz = (P_guess * self.nozzle.A_throat) / cstar
            return (m_dot_ox + m_dot_fuel) - mdot_noz

        # Solve for Pc
        sol = root_scalar(mass_balance, bracket=[1e5, 1e7], method="bisect")
        if not sol.converged:
            raise RuntimeError("Pc solve failed")

        P_cc = sol.root

        # Get thermo props at solved Pc
        T_c = self.C.get_Tcomb(P_cc, OF)
        MW, gamma = self.C.get_Chamber_MolWt_gamma(P_cc, OF, self.nozzle.expratio)
        R_spec = R_UNIV / MW

        # Nozzle outflow (you already have a solver for this)
        instThrust, _ = self.nozzle.sol_thrust(P_cc, T_c, gamma, R_spec)

        #setup for next timestep:
        self.r += r_dot*self.timestep

        print("out: ", P_cc, instThrust)

        return {"P_cc": P_cc, "thrust": instThrust}
