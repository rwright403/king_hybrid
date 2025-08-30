from src.models._tsw_n2o_constants import *
import numpy as np
from scipy.optimize import brentq


def secant(func, x1):
    x_eps = x1 * 0.005  # Set the tolerance to be 0.5% of init guess
    x2 = x1 -x1 * 0.01  # Set a second point 1% away from the original guess
    F1 = func(x1)  # Evaluate function at x1
    F2 = func(x2)  # Evaluate function at x2
    kk = 1  # Set up counter
    kk_max = 1000

    while np.abs(x2 - x1) >= (x_eps) and kk < kk_max:  # While error is too large and counter is less than max
        x3 = x2 - (F2 * (x2 - x1) / (F2 - F1)) 
        x1 = x2  # Move everything forward
        x2 = x3
        F1 = F2
        F2 = func(x2) 
        if (F1 == F2):
            return x2
        kk = kk + 1
    x = x2
    return x


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


def lightweight_span_wagner_eos_pressure(rho, T):
    tau = T_c / T
    delta = rho / rho_c
    ar_delta = np.sum(n1 * d1 * delta**(d1 - 1) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    return rho * R * T * (1 + delta * ar_delta)

def lightweight_span_wagner_eos_cp(rho, T):
    tau = T_c / T
    delta = rho / rho_c
    ao_tautau = -(c0 - 1) / tau**2 + np.sum(-v0 * u0**2 / T_c**2 * np.exp(-u0 * tau / T_c) / (1 - np.exp(-u0 * tau / T_c))**2)
    ar_tautau = np.sum(n1 * t1 * (t1 - 1) * tau**(t1 - 2) * delta**d1) + np.sum(n2 * t2 * (t2 - 2) * tau**(t2 - 2) * delta**d2 * np.exp(-delta**P0))
    ar_delta = np.sum(n1 * d1 * delta**(d1 - 1) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    ar_deltadelta = np.sum(n1 * d1 * (d1 - 1) * delta**(d1 - 2) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 2) * ((d2 - P0 * delta**P0) * (d2 - 1 - P0 * delta**P0) - P0**2 * delta**P0) * np.exp(-delta**P0))
    ar_deltatau = np.sum(n1 * d1 * t1 * delta**(d1 - 1) * tau**(t1 - 1)) + np.sum(n2 * t2 * tau**(t2 - 1) * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    
    num = (1 + delta*ar_delta - delta*tau*ar_deltatau)**2
    den = 1 + 2*delta*ar_delta + delta**2 * ar_deltadelta
    return R * (-tau**2 * (ao_tautau + ar_tautau) + num/den)

def lightweight_span_wagner_eos_d_rho_dT_P(rho, T):
    tau = T_c / T
    delta = rho / rho_c
    ar_delta = np.sum(n1 * d1 * delta**(d1 - 1) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    ar_deltadelta = np.sum(n1 * d1 * (d1 - 1) * delta**(d1 - 2) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 2) * ((d2 - P0 * delta**P0) * (d2 - 1 - P0 * delta**P0) - P0**2 * delta**P0) * np.exp(-delta**P0))
    ar_deltatau = np.sum(n1 * d1 * t1 * delta**(d1 - 1) * tau**(t1 - 1)) + np.sum(n2 * t2 * tau**(t2 - 1) * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))

    num = rho * R * (1 + delta*ar_delta - tau*delta*ar_deltatau)
    den = R*T * (1 + 2*delta*ar_delta + (delta**2)*ar_deltadelta)
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
        return self.rho * R * self.T * (1 + self.delta * self.ar_delta)

    # Specific internal energy [J/kg]
    @property
    def u(self):
        return R * self.T * self.tau * (self.ao_tau + self.ar_tau)

    # Specific entropy [J/kg/K]
    @property
    def s(self):
        return R * (self.tau * (self.ao_tau + self.ar_tau) - self.ao - self.ar)

    # Specific enthalpy [J/kg]
    @property
    def h(self):
        return R * self.T * (1 + self.tau*(self.ao_tau + self.ar_tau) + self.delta*self.ar_delta)

    # Helmholtz free energy per unit mass [J/kg]
    @property
    def a(self):
        return R * self.T * (self.ao + self.ar)

    # Chemical potential [J/kg]
    @property
    def mu(self):
        return self.a + self.p / self.rho

    # Isochoric specific heat [J/kg/K]
    @property
    def cv(self):
        return R * -self.tau**2 * (self.ao_tautau + self.ar_tautau)

    # Isobaric specific heat [J/kg/K]
    @property
    def cp(self):
        num = (1 + self.delta*self.ar_delta - self.delta*self.tau*self.ar_deltatau)**2
        den = 1 + 2*self.delta*self.ar_delta + self.delta**2 * self.ar_deltadelta
        return R * (-self.tau**2 * (self.ao_tautau + self.ar_tautau) + num/den)

    # du/drho|T
    @property
    def du_drho_const_T(self):
        return R * self.T / self.rho * (self.tau * self.delta * self.ar_deltatau)

    # dP/dT|rho
    @property
    def dP_dT_const_rho(self):
        return self.rho * R * (1 + self.delta*self.ar_delta - self.tau*self.delta*self.ar_deltatau)

    # dP/drho|T
    @property
    def dP_drho_const_T(self):
        return R * self.T * (1 + 2*self.delta*self.ar_delta + (self.delta**2)*self.ar_deltadelta)

    # d rho / dT |_P
    @property
    def d_rho_dT_const_P(self):
        num = self.rho * R * (1 + self.delta*self.ar_delta - self.tau*self.delta*self.ar_deltatau)
        den = R*self.T * (1 + 2*self.delta*self.ar_delta + (self.delta**2)*self.ar_deltadelta)
        return num / den
    


class SpanWagnerEOS_EquilibriumPhase(SpanWagnerEOS_BASE):
    """Equilibrium (saturation) properties using Span-Wagner EOS."""

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
            self.T = self.T_sat(P)

        # Find coexisting densities
        self.rho_v, self.rho_l = self.find_density_roots(self.T, self.P)

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
        return (P/(R*T))*(1/rho_v - 1/rho_l) - np.log(rho_l/rho_v) - (liq.ar - vap.ar)

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

    # -------------------------
    #   Public API
    # -------------------------
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
        return SpanWagnerEOS_SingleState(self.rho_l, self.T).h

    @property
    def h_sat_gas(self):
        return SpanWagnerEOS_SingleState(self.rho_v, self.T).h

    @property
    def rho_sat_liq(self):
        return self.rho_l

    @property
    def rho_sat_gas(self):
        return self.rho_v