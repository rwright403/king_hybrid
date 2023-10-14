#Adapted from Benjamin Klammer's MuleSim3 Matlab Script --> shoutout

#I love coolprop! ~ units: http://www.coolprop.org/v4/apidoc/CoolProp.html

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np

def secant(func, x1, U_tank, m_ox, V_tank):
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
        kk = kk + 1
    x = x2
    return x

def Verror(T_tank, U_tank, m_ox, V_tank):
    #tank temp evaluates to none when lambda used

    rho_liq = CP.PropsSI('D', 'Q', 0, 'T', T_tank, 'N2O')  # Density of liquid nitrous (kg/m^3)
    rho_vap = CP.PropsSI('D', 'Q', 1, 'T', T_tank, 'N2O')  # Density of nitrous gas (kg/m^3)

    u_liq = CP.PropsSI('U', 'Q', 0, 'T', T_tank, 'N2O')  # Internal Energy of liquid nitrous (J/kg)
    u_vap = CP.PropsSI('U', 'Q', 1, 'T', T_tank, 'N2O')  # Internal Energy of nitrous gas (J/kg)

    x_tank = (U_tank/m_ox - u_liq)/(u_vap - u_liq)
    V_tank_diff = m_ox * ( ((1-x_tank) /rho_liq) + (x_tank/rho_vap) ) - V_tank
    return V_tank_diff

def uerror(T, rho_tank, u_tank): #TODO:(u_tank + 7.3397e+5)
    u2 = thermo_span_wagner(rho_tank, T, 'u') + 7.3397e+5
    u_diff = u2 - u_tank
    return u_diff

def thermo_span_wagner(rho, T, param):
    # Constants for N2O
    R = 8.3144598 / 44.0128 * 1000  # Gas constant (J/kg*K)
    T_c = 309.52  # Critical Temperature (K)
    rho_c = 452.0115  # Critical Density (kg/m^3)

    n0 = np.array([0.88045, -2.4235, 0.38237, 0.068917, 0.00020367, 0.13122, 0.46032,
          -0.0036985, -0.23263, -0.00042859, -0.042810, -0.023038])
    n1 = n0[0:5]
    n2 = n0[5:12]
    a1 = 10.7927224829
    a2 = -8.2418318753
    c0 = 3.5
    v0 = np.array([2.1769, 1.6145, 0.48393])
    u0 = np.array([879, 2372, 5447])
    t0 = np.array([0.25, 1.125, 1.5, 0.25, 0.875, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5])
    d0 = np.array([1, 1, 1, 3, 7, 1, 2, 5, 1, 1, 4, 2])
    P0 = np.array([1, 1, 1, 2, 2, 2, 3])
    t1 = t0[0:5]
    t2 = t0[5:12]
    d1 = d0[0:5]
    d2 = d0[5:12]

    # Calculate non-dimensional variables
    tau = T_c / T
    delta = rho / rho_c

    # Calculate explicit Helmholtz energy and derivatives
    ao = a1 + a2 * tau + np.log(delta) + (c0 - 1) * np.log(tau) + np.sum(v0 * np.log(1 - np.exp(-u0 * tau / T_c)))
    ar = np.sum(n1 * tau**t1 * delta**d1) + np.sum(n2 * tau**t2 * delta**d2 * np.exp(-delta**P0))
    ao_tau = a2 + (c0 - 1) / tau + np.sum(v0 * u0 / T_c * np.exp(-u0 * tau / T_c) / (1 - np.exp(-u0 * tau / T_c)))
    ao_tautau = -(c0 - 1) / tau**2 + np.sum(-v0 * u0**2 / T_c**2 * np.exp(-u0 * tau / T_c) / (1 - np.exp(-u0 * tau / T_c))**2)
    ar_tau = np.sum(n1 * t1 * tau**(t1 - 1) * delta**d1) + np.sum(n2 * t2 * tau**(t2 - 1) * delta**d2 * np.exp(-delta**P0))
    ar_tautau = np.sum(n1 * t1 * (t1 - 1) * tau**(t1 - 2) * delta**d1) + np.sum(n2 * t2 * (t2 - 2) * tau**(t2 - 2) * delta**d2 * np.exp(-delta**P0))
    ar_delta = np.sum(n1 * d1 * delta**(d1 - 1) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    ar_deltadelta = np.sum(n1 * d1 * (d1 - 1) * delta**(d1 - 2) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 2) * ((d2 - P0 * delta**P0) * (d2 - 1 - P0 * delta**P0) - P0**2 * delta**P0) * np.exp(-delta**P0))
    ar_deltatau = np.sum(n1 * d1 * t1 * delta**(d1 - 1) * tau**(t1 - 1)) + np.sum(n2 * t2 * tau**(t2 - 1) * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))

    out = 0.0
    if param == 'p':  # Pressure (Pa)
        out = rho * R * T * (1 + delta * ar_delta)
    elif param == 'u':  # Specific internal energy (J/kg)
        out = R * T * tau * (ao_tau + ar_tau)
    elif param == 's':  # Specific entropy (J/kg*K)
        out = R * (tau * (ao_tau + ar_tau) - ao - ar)
    elif param == 'h':  # Specific enthalpy (J/kg)
        out = R * T * (1 + tau * (ao_tau + ar_tau) + delta * ar_delta)
    elif param == 'cv':  # Specific heat constant volume (J/kg*K)
        out = R * -tau**2 * (ao_tautau + ar_tautau)
    elif param == 'cp':  # Specific heat constant pressure (J/kg*K)
        out = R * (-tau**2 * (ao_tautau + ar_tautau) + (1 + delta * ar_delta - delta * tau * ar_deltatau)**2 / (1 + 2 * delta * ar_delta + delta**2 * ar_deltadelta))
    elif param == 'a':  # Speed of sound (m/s)
        out = np.sqrt(R * T * (1 + 2 * delta * ar_delta + delta**2 * ar_deltadelta - (1 + delta * ar_delta - delta * tau * ar_deltatau)**2 / (tau**2 * (ao_tautau + ar_tautau))))
    else:
        raise ValueError('Invalid input')

    return out

class OxTank():
    def __init__(self, oxidizer, timestep, fill_level, C_inj, V_tank, P_tank, P_cc, all_error):
        self.oxidizer = oxidizer
        self.timestep = timestep
        self.fill_level = fill_level
        self.C_inj = C_inj
        self.m_dot_ox = 0
        self.V_tank = V_tank
        self.P_tank = P_tank
        self.P_cc = P_cc
        self.all_error = all_error

        #setup 
        self.t = 0
        self.m_dot_ox_prev = self.m_dot_ox

        #setup thermo data
        self.rho_liq = CP.PropsSI('D', 'Q', 0, 'P', self.P_tank, 'N2O')  # Density of liquid nitrous (kg/m^3)
        self.rho_vap = CP.PropsSI('D', 'Q', 1, 'P', self.P_tank, 'N2O')  # Density of nitrous gas (kg/m^3)

        self.u_liq = CP.PropsSI('U', 'Q', 0, 'P', self.P_tank, 'N2O')  # Internal Energy of liquid nitrous (kg/m^3)
        self.u_vap = CP.PropsSI('U', 'Q', 1, 'P', self.P_tank, 'N2O')  # Internal Energy of nitrous gas (kg/m^3)

        #Calculate m_ox from fill level
        self.m_ox =  self.V_tank*(self.rho_liq*self.fill_level+self.rho_vap*(1-self.fill_level))

        self.x_tank = ( (self.V_tank/self.m_ox) - ((self.rho_liq)**-1) )/( ((self.rho_vap)**-1) - ((self.rho_liq)**-1)) #quality

        self.T_tank = CP.PropsSI('T', 'Q', self.x_tank, 'P', self.P_tank, 'N2O')  # Temperature of nitrous gas (kg/m^3)

        self.u_tank = self.x_tank*self.u_vap + (1-self.x_tank)*self.u_liq
        self.U_tank = self.m_ox*self.u_tank
          

    def inst(self,P_cc):
        self.P_cc = P_cc

        #setup iteration error tolerance
        self.V_tank_err = self.all_error*self.V_tank
        self.u_tank_err = self.all_error*self.u_tank

        ### start
        if self.x_tank < 1:
            while np.abs(Verror(self.T_tank, self.U_tank, self.m_ox, self.V_tank ) ) > self.V_tank_err:
                self.T_tank = secant((lambda T: Verror(T, self.U_tank, self.m_ox, self.V_tank)), self.T_tank, self.U_tank, self.m_ox, self.V_tank)

            #use temperature to calculate thermo properties of tank
            self.P_tank = CP.PropsSI('P', 'Q', self.x_tank, 'T', self.T_tank, 'N2O')
            h_liq = CP.PropsSI('H', 'Q', 0, 'T', self.T_tank, 'N2O')
            h_vap = CP.PropsSI('H', 'Q', 1, 'T', self.T_tank, 'N2O')
            self.rho_liq = CP.PropsSI('D', 'Q', 0, 'T', self.T_tank, 'N2O')
            self.rho_vap = CP.PropsSI('D', 'Q', 1, 'T', self.T_tank, 'N2O')
            self.u_liq = CP.PropsSI('U', 'Q', 0, 'T', self.T_tank, 'N2O')
            self.u_vap = CP.PropsSI('U', 'Q', 1, 'T', self.T_tank, 'N2O')

            self.x_tank = (self.U_tank/self.m_ox - self.u_liq)/(self.u_vap - self.u_liq)
            self.u_tank = self.x_tank*self.u_vap + (1 - self.x_tank)*self.u_liq
            h_tank = self.x_tank*h_vap + (1 - self.x_tank)*h_liq
            self.rho_tank = self.x_tank*self.rho_vap + (1-self.x_tank)*self.rho_liq

            #update current time
            self.t = self.t + self.timestep
            #assume only liquid draining from tank
            self.rho_exit = self.rho_liq
            h_exit = h_liq

        else:
            #solve variables
            self.rho_tank = self.m_ox/self.V_tank
            self.u_tank = self.U_tank/self.m_ox

            #TODO: uerror frunction
            while np.abs(uerror(self.T_tank, self.rho_tank, self.u_tank) ) > self.u_tank_err:
                self.T_tank = secant((lambda T: uerror(T, self.rho_tank, self.u_tank)), self.T_tank, self.U_tank, self.m_ox, self.V_tank)

            self.P_tank = thermo_span_wagner(self.rho_tank, self.T_tank, 'p')
            h_tank = thermo_span_wagner(self.rho_tank, self.T_tank, 'u')

            #update current time
            self.t = self.t + self.timestep
            #rho exit is rho vapor (assume only vapor left in tank)
            self.rho_exit = self.rho_tank
            h_exit = h_tank + 7.3397e+05 #Convert from Span-Wagner enthalpy convention to NIST


        #injector model
        self.m_dot_ox = self.C_inj*np.sqrt(2*self.rho_exit*(self.P_tank -self.P_cc)) #this uses a incompressible fluid assumption #TODO: add feed pressure
        #TODO: BELOW IS A PLACEHOLDER: INTEGRATE WITH THRUST CURVE
        #self.P_cc = self.P_cc*0.9 #Pa

        #Ben does this to eliminate numerical instability
        if self.t == self.timestep:
            self.m_dot_ox = 0.5 * self.m_dot_ox
        else:
            self.m_dot_ox = 0.5 * self.m_dot_ox + 0.5 * self.m_dot_ox_prev

        #move forward in time with differential eqns
        self.m_ox = self.m_ox - self.m_dot_ox*self.timestep
        self.m_dot_ox_prev = self.m_dot_ox
        self.U_tank = self.U_tank -self.m_dot_ox*h_exit*self.timestep