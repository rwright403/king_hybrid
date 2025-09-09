import numpy as np
from scipy.optimize import brentq
from src.utils.numerical_methods import secant
from src.models._thermo.n2o_thermo_span_wagner_n2o_constants import *

REF_SHIFT = 7.3397e+05 #Convert from Span-Wagner enthalpy convention to NIST


import CoolProp.CoolProp as CP


class SpanWagnerEOS_BASE:

    def compute_helmholtz_energy_and_derivs(self): 
        # Calculate explicit Helmholtz energy and derivatives
        self.ao = a1 + a2 * self.tau + np.log(self.delta) + (c0 - 1) * np.log(self.tau) + np.sum(v0 * np.log(1 - np.exp(-u0 * self.tau / T_c)))
        self.ar = np.sum(n1 * self.tau**t1 * self.delta**d1) + np.sum(n2 * self.tau**t2 * self.delta**d2 * np.exp(-self.delta**P0))
        self.ao_tau = a2 + (c0 - 1) / self.tau + np.sum(v0 * u0 / T_c * np.exp(-u0 * self.tau / T_c) / (1 - np.exp(-u0 * self.tau / T_c)))
        self.ao_tautau = -(c0 - 1) / self.tau**2 + np.sum(-v0 * u0**2 / T_c**2 * np.exp(-u0 * self.tau / T_c) / (1 - np.exp(-u0 * self.tau / T_c))**2)
        self.ar_tau = np.sum(n1 * t1 * self.tau**(t1 - 1) * self.delta**d1) + np.sum(n2 * t2 * self.tau**(t2 - 1) * self.delta**d2 * np.exp(-self.delta**P0))
        self.ar_tautau = np.sum(n1 * t1 * (t1 - 1) * self.tau**(t1 - 2) * self.delta**d1) + np.sum(n2 * t2 * (t2 - 2) * self.tau**(t2 - 2) * self.delta**d2 * np.exp(-self.delta**P0))
        self.ar_delta = np.sum(n1 * d1 * self.delta**(d1 - 1) * self.tau**t1) + np.sum(n2 * self.tau**t2 * self.delta**(d2 - 1) * (d2 - P0 * self.delta**P0) * np.exp(-self.delta**P0))
        self.ar_deltadelta = np.sum(n1 * d1 * (d1 - 1) * self.delta**(d1 - 2) * self.tau**t1) + np.sum(n2 * self.tau**t2 * self.delta**(d2 - 2) * ((d2 - P0 * self.delta**P0) * (d2 - 1 - P0 * self.delta**P0) - P0**2 * self.delta**P0) * np.exp(-self.delta**P0))
        self.ar_deltatau = np.sum(n1 * d1 * t1 * self.delta**(d1 - 1) * self.tau**(t1 - 1)) + np.sum(n2 * t2 * self.tau**(t2 - 1) * self.delta**(d2 - 1) * (d2 - P0 * self.delta**P0) * np.exp(-self.delta**P0))


def lightweight_span_wagner_eos_pressure(rho, T): #NOTE: RETURNS P IN [Pa]
    tau = T_c / T
    delta = rho / rho_c
    ar_delta = np.sum(n1 * d1 * delta**(d1 - 1) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    return rho * R_mass * T * (1 + delta * ar_delta)

def lightweight_span_wagner_eos_cp(rho, T): #NOTE: RETURNS CP IN [J/(kg K)]
    tau = T_c / T
    delta = rho / rho_c
    ao_tautau = -(c0 - 1) / tau**2 + np.sum(-v0 * u0**2 / T_c**2 * np.exp(-u0 * tau / T_c) / (1 - np.exp(-u0 * tau / T_c))**2)
    ar_tautau = np.sum(n1 * t1 * (t1 - 1) * tau**(t1 - 2) * delta**d1) + np.sum(n2 * t2 * (t2 - 2) * tau**(t2 - 2) * delta**d2 * np.exp(-delta**P0))
    ar_delta = np.sum(n1 * d1 * delta**(d1 - 1) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    ar_deltadelta = np.sum(n1 * d1 * (d1 - 1) * delta**(d1 - 2) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 2) * ((d2 - P0 * delta**P0) * (d2 - 1 - P0 * delta**P0) - P0**2 * delta**P0) * np.exp(-delta**P0))
    ar_deltatau = np.sum(n1 * d1 * t1 * delta**(d1 - 1) * tau**(t1 - 1)) + np.sum(n2 * t2 * tau**(t2 - 1) * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    
    num = (1 + delta*ar_delta - delta*tau*ar_deltatau)**2
    den = 1 + 2*delta*ar_delta + delta**2 * ar_deltadelta
    #NOTE: RETURNS CP IN [J/(kg K)]
    return R_mass * (-tau**2 * (ao_tautau + ar_tautau) + num/den) 

def lightweight_span_wagner_eos_d_rho_dT_P(rho, T): #NOTE: RETURNS D_RHO_DT_P IN [(kg/m^3)/K]
    tau = T_c / T
    delta = rho / rho_c
    ar_delta = np.sum(n1 * d1 * delta**(d1 - 1) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    ar_deltadelta = np.sum(n1 * d1 * (d1 - 1) * delta**(d1 - 2) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 2) * ((d2 - P0 * delta**P0) * (d2 - 1 - P0 * delta**P0) - P0**2 * delta**P0) * np.exp(-delta**P0))
    ar_deltatau = np.sum(n1 * d1 * t1 * delta**(d1 - 1) * tau**(t1 - 1)) + np.sum(n2 * t2 * tau**(t2 - 1) * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))

    num = rho * R_mass * (1 + delta*ar_delta - tau*delta*ar_deltatau)
    den = R_mass * T * (1 + 2*delta*ar_delta + (delta**2)*ar_deltadelta)
    return num / den 


class SpanWagnerEOS_SingleState(SpanWagnerEOS_BASE):

    def __init__(self, rho, T):
        self.rho = rho
        self.T = T
        self.tau = T_c / T
        self.delta = rho / rho_c
        self.compute_helmholtz_energy_and_derivs() 

    # Pressure [Pa]
    @property
    def P(self):
        return self.rho * R_mass * self.T * (1 + self.delta * self.ar_delta)

    # Specific internal energy [J/kg]
    @property
    def u(self):
        return R_mass * self.T * self.tau * (self.ao_tau + self.ar_tau) + REF_SHIFT

    # Specific enthalpy [J/kg]
    @property
    def h(self):
        return R_mass * self.T * (1 + self.tau*(self.ao_tau + self.ar_tau) + self.delta*self.ar_delta) + REF_SHIFT

    # Specific entropy [J/kg/K]
    @property
    def s(self):
        return R_mass * (self.tau * (self.ao_tau + self.ar_tau) - self.ao - self.ar)

    # Helmholtz free energy per unit mass [J/kg]
    @property
    def a(self):
        return R_mass * self.T * (self.ao + self.ar)

    # Chemical potential [J/kg]
    @property
    def mu(self):
        return self.a + self.P / self.rho

    # Isochoric specific heat [J/kg/K]
    @property
    def cv(self):
        return R_mass * -self.tau**2 * (self.ao_tautau + self.ar_tautau)

    # Isobaric specific heat [J/kg/K]
    @property
    def cp(self):
        num = (1 + self.delta*self.ar_delta - self.delta*self.tau*self.ar_deltatau)**2
        den = 1 + 2*self.delta*self.ar_delta + self.delta**2 * self.ar_deltadelta
        return R_mass * (-self.tau**2 * (self.ao_tautau + self.ar_tautau) + num/den)

    # du/drho|T [(J*m^3)/kg^2]
    @property
    def du_drho_const_T(self):
        return R_mass * self.T / self.rho * (self.tau * self.delta * self.ar_deltatau)

    # dP/dT|rho [kg/(m^3*K)]
    @property
    def dP_dT_const_rho(self):
        return self.rho * R_mass * (1 + self.delta*self.ar_delta - self.tau*self.delta*self.ar_deltatau)

    # dP/drho|T [(Pa*m^3)/kg]
    @property
    def dP_drho_const_T(self):
        return R_mass * self.T * (1 + 2*self.delta*self.ar_delta + (self.delta**2)*self.ar_deltadelta)

    # drho/dT|_P [kg/(m^3*K)]
    @property
    def d_rho_dT_const_P(self):
        num = self.rho * R_mass * (1 + self.delta*self.ar_delta - self.tau*self.delta*self.ar_deltatau)
        den = R_mass * self.T * (1 + 2*self.delta*self.ar_delta + (self.delta**2)*self.ar_deltadelta)
        return num / den
    


class SpanWagnerEOS_EquilibriumPhase(SpanWagnerEOS_BASE):
    """NOTE: NOT FULLY THERMODYNAMICALLY CLOSED, look here for where I cheat and use COOLPROP. Equilibrium (saturation) properties using Span-Wagner EOS."""

    def __init__(self, T=None, P=None):
        if T is None and P is None:
            raise ValueError("Must specify either T or P")

        self.T = T
        self.P = P

        # If T given, solve P_sat
        if T is not None and P is None:
            self.P = self.P_sat(T)

        # If P given, solve T_sat
        if P is not None and T is None:
#NOTE: TESTING W COOLPROP:
            # OLD self.T = self.T_sat(P)
            self.T = CP.PropsSI('T', 'P', P, 'Q', 1, "N2O")


        # Find coexisting densities
        self.rho_v, self.rho_l = self.find_density_roots(self.T, self.P) #NOTE: UNITS [kg/m^3]

    # -------------------------
    #   Core saturation solvers
    # -------------------------
    @staticmethod
    def P_sat_anc(T):
        A, B, C = 4.80716087, 967.819748, 19.6368887
        if 140 < T < 310:
            return 1e5 * (10**(A - (B / (T+C))))
        raise ValueError("T outside range")

    @staticmethod
    def T_sat_anc(P):
        A, B, C = 4.80716087, 967.819748, 19.6368887
        T_sat_est = B/(A - np.log10(P/1e5)) - C
        if 140 < T_sat_est < 310:
            return T_sat_est
        raise ValueError("T outside range")


    @staticmethod
    def span2000_residual(P, T):
        rho_v, rho_l = SpanWagnerEOS_EquilibriumPhase.find_density_roots(T, P)
        vap = SpanWagnerEOS_SingleState(rho_v, T)
        liq = SpanWagnerEOS_SingleState(rho_l, T)
        return (P/(R_mass*T))*(1/rho_v - 1/rho_l) - np.log(rho_l/rho_v) - (liq.ar - vap.ar)

    @staticmethod
    def find_density_roots(T, P_target, rho_min=0.1, rho_max=1200, n_points=50):
        rhos = np.linspace(rho_min, rho_max, n_points)
        residuals = [SpanWagnerEOS_SingleState(r,T).P - P_target for r in rhos]
        signs = np.sign(residuals)
        crossings = np.where(np.diff(signs))[0]

        roots = []
        for i in crossings:
            try:
                root = brentq(
                    lambda rho: SpanWagnerEOS_SingleState(rho,T).P - P_target,
                    rhos[i], rhos[i+1]
                )
                roots.append(root)
            except ValueError:
                pass

        if len(roots) < 2:
            raise RuntimeError("Less than two density roots found")
        return sorted(roots)[0], sorted(roots)[-1]  # (vap, liq)


    @classmethod
    def P_sat(cls, T):
        P_est = SpanWagnerEOS_EquilibriumPhase.P_sat_anc(T)  # ancillary guess
        pcnt = 0.1
        while pcnt < 0.5:
            try:
                P_min, P_max = (1-pcnt)*P_est, (1+pcnt)*P_est
                return brentq(lambda P: cls.span2000_residual(P, T),
                              P_min, P_max, xtol=1e-6, rtol=1e-6)
            except ValueError:
                pcnt += 0.05
        raise RuntimeError("P_sat solver failed")

    @classmethod
    def T_sat(cls, P):
        T_guess = SpanWagnerEOS_EquilibriumPhase.T_sat_anc(P)
        return secant(lambda T: cls.P_sat(T) - P, T_guess)

    @property
    def h_sat_liq(self):
        return SpanWagnerEOS_SingleState(self.rho_l, self.T).h # already added REF_SHIFT in this call

    @property
    def h_sat_gas(self):
        return SpanWagnerEOS_SingleState(self.rho_v, self.T).h # already added REF_SHIFT in this call

    @property
    def u_sat_liq(self):
        return SpanWagnerEOS_SingleState(self.rho_l, self.T).u # already added REF_SHIFT in this call

    @property
    def u_sat_gas(self):
        return SpanWagnerEOS_SingleState(self.rho_v, self.T).u # already added REF_SHIFT in this call

    @property
    def rho_sat_liq(self):
        return self.rho_l

    @property
    def rho_sat_gas(self):
        return self.rho_v
    


### Find Spinodal:

def lightweight_dP_drho_constT(rho, T):
    
    R_u = 8.314462618 # [J/(mol-K)]
    MM = 0.0440128 # [kg/mol]

    R_mass = R_u/MM # Specific Gas constant [J/(kg*K)]


    tau = T_c / T
    delta = rho / rho_c
    ar_delta = np.sum(n1 * d1 * delta**(d1 - 1) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    ar_deltadelta = np.sum(n1 * d1 * (d1 - 1) * delta**(d1 - 2) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 2) * ((d2 - P0 * delta**P0) * (d2 - 1 - P0 * delta**P0) - P0**2 * delta**P0) * np.exp(-delta**P0))

    return R_mass * T * (1 + 2*delta*ar_delta + (delta**2)*ar_deltadelta)

def find_spinodal(T, rho_c=452.0115 , rho_min=1e-3, rho_max=2000.0):
    
    # vapor-side spinodal
    rho_v_sp = brentq(
        lambda r: lightweight_dP_drho_constT(r, T),
        rho_min, 0.9 * rho_c
    )
    # liquid-side spinodal
    rho_l_sp = brentq(
        lambda r: lightweight_dP_drho_constT(r, T),
        1.1 * rho_c, rho_max
    )
    return rho_v_sp, rho_l_sp

