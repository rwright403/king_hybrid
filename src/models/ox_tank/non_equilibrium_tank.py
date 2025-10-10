import numpy as np
from scipy.optimize import root, root_scalar, brentq
from src.models.ox_tank._base import BaseTank
from thermo import Chemical
import matplotlib.pyplot as plt
import traceback

from src.models.inj._base import build_state

from src.models._thermo.convection_heat_transfer import *
from src.models._thermo.conduction_heat_transfer import *

from src.models._thermo.n2o_thermo_span_wagner_class import SpanWagnerEOS_SingleState, SpanWagnerEOS_EquilibriumPhase, lightweight_span_wagner_eos_pressure, lightweight_span_wagner_eos_cp, lightweight_span_wagner_eos_d_rho_dT_P

from src.utils.numerical_methods import rk4_step


# Global Constants:
R_U = 8.31446 #J/(mol K)
T_REF = 298.15 #K

n2o_g = Chemical('nitrous oxide', T=T_REF)
PC = n2o_g.Pc
TC = n2o_g.Tc
OMEGA = n2o_g.omega

MW = (n2o_g.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o_g.omega - 0.26992*n2o_g.omega**2
b = 0.07780*(R_U*TC/PC)
g = 9.81 #m/s^2

CW = 896 #J/(kg K) ~ this is the specific heat capacity of the wall material, for Al 6061 from MATWEB: https://www.matweb.com/search/datasheet.aspx?MatGUID=b8d536e0b9b54bd7b69e4124d8f1d20a&ckck=1

E = 2.1e4 #1.3e3 #2.1e4 #NOTE: EMPIRICAL FACTOR E for high mass flow effects to scale Q_liq_to_sat_surf and m_dot_evap. This was obtained from [7], which claims E is an oxidizer property. This is contradicted by [8] which claims E is a function of the tank as a system


"""
#TODO: update with other injector model once we get this thing up
def spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit):


    #for tomacz test case!!!
"""
"""
    pipe_inj_time = [ 0, 0.25, 1.1, 1.5,4]
    pipe_inj_m_dot = [ (-50/1000), (-43/1000), (-41.8/1000), (-36/1000), (-22/1000)]

    m_dot_spi = np.interp(LOOKUP_TIME, pipe_inj_time , pipe_inj_m_dot)
"""
"""
    #NOTE: PIPING M_DOT TO ISOLATE REST OF MODEL FOR DEBUG
    m_dot_spi = -3.75
    return m_dot_spi
"""

def solve_m_dot_evap(liq_state, sat_surf, Q_dot_liq_to_sat_surf, Q_dot_sat_surf_to_gas):
    m_dot_evap = 0
    if (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) > 0:

        m_dot_evap = (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) / ( (sat_surf.h_sat_gas-sat_surf.h_sat_liq) + (sat_surf.h_sat_liq - liq_state.h)  ) 

    return m_dot_evap


"""
def solve_m_dot_condensed(sat_surf, gas_state, V_gas):
    m_dot_cond = 0

    if (P_tank > sat_surf.P):
        m_dot_cond = ((gas_state.P-sat_surf.P)*V_gas*MW)/( (R_U/MW)*gas_state.T*(TIMESTEP) )   #NOTE EDIT DENOM FOR TESTING, OLD FOR REF: ( preos_g.Z_g*(R_U/MW)*T_gas*(TIMESTEP) )

    #NOTE: CONDENSATION FROM PREOS
    #if p_tank > p_sat_gas, then condensation to enforce equilibrium

    ###NOTE: try m_dot_cond = 0 FOR DEBUGGING #NOTE
    m_dot_cond = 0
    return m_dot_cond
"""


def solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj):
    m_dot_liq = -m_dot_evap + m_dot_cond + m_dot_inj
    m_dot_gas = m_dot_evap - m_dot_cond #convert sign convention from liq cv to gas cv

    return m_dot_liq, m_dot_gas


def thermo_residuals(rhos, T_liq, T_gas, m_liq, m_gas, V_tank):
    rho_liq, rho_gas = rhos

    P_liq = lightweight_span_wagner_eos_pressure(rho_liq, T_liq)
    P_gas = lightweight_span_wagner_eos_pressure(rho_gas, T_gas)

    V_est = (m_liq / rho_liq) + (m_gas / rho_gas)

    #print("P liq, gas, V_est, diff: ", P_liq, P_gas, V_est, (P_liq - P_gas) )

    return [
        P_liq - P_gas,      # pressure equilibrium
        V_est - V_tank      # volume constraint
    ]

def solve_thermo_params(T_liq, T_gas, m_liq, m_gas, rho_liq_prev, rho_gas_prev, V_tank):
    # Initial guess for densities: previous timestep values or saturation values
    guess = [rho_liq_prev, rho_gas_prev]

    try:
        sol = root(thermo_residuals, guess, args=(T_liq, T_gas, m_liq, m_gas, V_tank), method='hybr')
    except ValueError as e:
        print(e)

    """if not sol.success or np.linalg.norm(sol.fun) > volume_err_tolerance:
        raise RuntimeError("solve_thermo_params: Convergence failed")"""

    rho_liq, rho_gas = sol.x

    # Calculate common pressure
    P_tank = lightweight_span_wagner_eos_pressure(rho_gas, T_gas)

    return rho_liq, rho_gas, P_tank




def single_solve_T_dot_liq_gas(V_dot_liq, liq_state, sat_surf, gas_state, P_tank, m_liq, m_gas, V_liq, V_gas, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, debug_mode):

    m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj)

    V_dot_gas = -V_dot_liq

    d_rho_dt_liq = (1/V_liq)*m_dot_liq -(m_liq/(V_liq**2))*V_dot_liq
    d_rho_dt_gas = (1/V_gas)*m_dot_gas -(m_gas/(V_gas**2))*V_dot_gas


    # NOTE: BROKE FOR TESTING ON PURPOSE  no cond!!!!!

    # abs
    U_dot_liq = m_dot_inj*liq_state.h - m_dot_evap*liq_state.h + m_dot_cond*(sat_surf.h_sat_liq) - P_tank*V_dot_liq + Q_dot_liq
    U_dot_gas =                         m_dot_evap*(sat_surf.h_sat_gas) - m_dot_cond*(gas_state.h) - P_tank*V_dot_gas + Q_dot_gas


    T_dot_liq = (1/liq_state.cv)*( (1/m_liq) * (U_dot_liq - (liq_state.u * m_dot_liq)) - (liq_state.du_drho_const_T * d_rho_dt_liq) )
    T_dot_gas = (1/gas_state.cv)*( (1/m_gas) * (U_dot_gas - (gas_state.u * m_dot_gas)) - (gas_state.du_drho_const_T * d_rho_dt_gas) )


    if debug_mode == True:
        a=1
        #print("T_dot_gas: ", T_dot_gas, U_dot_gas,  m_dot_evap*(sat_surf.h_sat_gas), - m_dot_cond*(gas_state.h), - P_tank*V_dot_gas, + Q_dot_gas)
        #print("T_dot_gas: ", T_dot_gas, (1/gas_state.cv)*(1/m_gas)*(U_dot_gas - (gas_state.u * m_dot_gas)) , -(1/gas_state.cv)*(gas_state.du_drho_const_T * d_rho_dt_gas) )


    return T_dot_liq, T_dot_gas 

def P_dot_error(V_dot_guess, liq_state, sat_surf, gas_state, P_tank, m_liq, m_gas, V_liq, V_gas, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas):    

    V_dot_gas = -V_dot_guess #guessing for liquid

    m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj) #or put this outside and pass in?


    d_rho_dt_liq = (1/V_liq)*m_dot_liq - (m_liq/(V_liq**2))*V_dot_guess
    d_rho_dt_gas = (1/V_gas)*m_dot_gas - (m_gas/(V_gas**2))*V_dot_gas

    T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_guess, liq_state, sat_surf, gas_state, P_tank, m_liq, m_gas, V_liq, V_gas, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, False)

    P_dot_liq = liq_state.dP_dT_const_rho*T_dot_liq + liq_state.dP_drho_const_T*d_rho_dt_liq

    P_dot_gas = gas_state.dP_dT_const_rho*T_dot_gas + gas_state.dP_drho_const_T*d_rho_dt_gas

    return P_dot_liq - P_dot_gas


class non_equilibrium_tank_model(BaseTank):
    def __init__(self, timestep, m_ox, P_tank, P_atm, T_atm, rho_atm, V_tank, diam_out, diam_in, rho_wall, k_w, volume_err_tol, P_dot_err_tol, injector):
        super().__init__(injector, timestep)
        self.TIMESTEP = timestep
        self.volume_err_tolerance = volume_err_tol
        self.P_dot_err_tolerance = P_dot_err_tol

        self.timestep = timestep

        self.T_atm = T_atm
        self.P_atm = P_atm
        self.rho_atm = rho_atm

        #NOTE: Tank geometry and setup
        self.rho_wall = rho_wall
        self.k_w = k_w
        self.diam_out = diam_out
        self.diam_in = diam_in
        self.Tank_Inner_Area = 0.25*np.pi*(self.diam_in**2)
        self.Tank_xsec_Area = 0.25*np.pi*((self.diam_out**2)-(self.diam_in**2))

        if self.diam_in >= self.diam_out:
            raise ValueError("Tank inner diameter >= Tank outer diameter, this will give nonphysical results so check your inputs!")

        #NOTE: System starts with valve closed assuming thermal equillibrium:
        self.m_dot_liq = 0
        self.m_dot_gas = 0

        self.P_cc = self.P_atm

        self.P_tank = P_tank

        #### Setup Thermo Properties - Assuming Tank Starts at thermal equilibrium --> Sat Conditions

        self.V_tank = V_tank
        self.m_ox = m_ox

        initial_eq_state = SpanWagnerEOS_EquilibriumPhase(None, P_tank)

        self.rho_liq = initial_eq_state.rho_sat_liq
        self.rho_gas = initial_eq_state.rho_sat_gas
        
        rho_bulk_tank = self.m_ox/V_tank
        x_tank = ( (1/rho_bulk_tank)-(1/self.rho_liq) ) / ( (1/self.rho_gas)-(1/self.rho_liq) )


        #gas cv setup
        self.T_gas = T_atm
        self.m_gas = x_tank*self.m_ox
        v_gas = 1/self.rho_gas
        self.V_gas = v_gas*self.m_gas

        #liquid cv
        self.T_liq = T_atm
        self.m_liq = self.m_ox-self.m_gas
        v_liq = 1/self.rho_liq
        self.V_liq = v_liq*self.m_liq


        self.V_tank = self.V_liq+self.V_gas # resolving just to get rid of unstability "what are you going to do if the aluminum is too small? water it? give it sunlight? let it grow?"
        self.height_tank = self.V_tank/(0.25*np.pi*(diam_in**2))

        # NOTE: ASSUME TANK WALL AT T_ATM, MAKE SURE THIS MATCHES REAL WORLD INITIAL CONDITIONS
        self.T_wall_gas = T_atm 
        self.T_wall_liq = T_atm


        self.rho_liq_prev = self.rho_liq
        self.rho_gas_prev = self.rho_gas

        self.V_dot_liq_prev = 1E-8

        self.t = 0.0
        self.m_dot_ox_prev = 0.0


        self.m_tank = 0.25*np.pi*(self.diam_out**2-self.diam_in**2) * self.height_tank

        print("init non eq tank: ", x_tank, self.V_tank, self.rho_liq, rho_bulk_tank, self.rho_gas)

        self.state = build_state()
   



    def system_of_liq_odes(self, t, y, P_cc):

        T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas = y  # Unpack state variables

        ### Solve thermo parameters! - old according to [8]
        rho_liq, rho_gas, P_tank = solve_thermo_params(T_liq, T_gas, m_liq, m_gas, self.rho_liq_prev, self.rho_gas_prev, self.V_tank)

        gas_state = SpanWagnerEOS_SingleState(rho_gas, T_gas) # I get this might be clunky but I was having serious speed issues so I tried a lightweight pressure sol in the iterative method. Yes, there are better ways to setup this entire program in general, my bad, im a lot older and more pythonic than when i started this script
        sat_surf = SpanWagnerEOS_EquilibriumPhase(None, P_tank)
        liq_state = SpanWagnerEOS_SingleState(rho_liq, T_liq)
        
        # Mass transfer (1) from injector
        
        self.state["P_1"] = P_tank
        self.state["P_sat"] = P_tank
        self.state["P_2"] = P_cc
        self.state["rho_1"] = rho_liq
        self.state["h_1"] = liq_state.h
        self.state["T_1"] = T_liq

        m_dot_inj = -self.injector.m_dot(self.state)
        # Apply smoothing like equilibrium tank
        if self.t <= (self.timestep*2):   # first step
            m_dot_inj = 0.5 * m_dot_inj
        else:
            m_dot_inj = 0.5 * m_dot_inj + 0.5 * (-self.m_dot_ox_prev)


        # Heat transfer (2) from saturated surface to gas                       (T_1, T_2, P_tank, rho_2, c, n, tank_diam, fluid)
        # L = tank inner diam , Area of circle x section
        T_film_gas = ((sat_surf.T + T_gas)/2 )
        Q_dot_sat_surf_to_gas = solve_Q_dot_natural_convection_gas(rho_gas, sat_surf.T, T_gas, T_film_gas, P_tank, 0.15, 0.333, self.diam_in, self.Tank_Inner_Area, "N2O" ) #relative to gas cv
    
        # Heat transfer (3)  from liq to saturated surface (sat surface assumed to be a liquid with quality 0)
        T_film_liq = ((sat_surf.T + T_liq)/2 )
        Q_dot_liq_to_sat_surf = (E)*solve_Q_dot_natural_convection_liq(rho_liq, T_liq, sat_surf.T, T_film_liq, P_tank, 0.15, 0.333, self.diam_in, self.Tank_Inner_Area, "N2O" ) #relative to liq cv
        #NOTE:CORRECTION FACTOR for nitrous oxide heat transfer E = (E) to account for blowing as per [7],[8]

        # Mass transfer (3) by evaporation 
        m_dot_evap = solve_m_dot_evap(liq_state, sat_surf, Q_dot_liq_to_sat_surf, Q_dot_sat_surf_to_gas)

        # Mass transfer (2) by condensation
        V_gas = m_gas/rho_gas
        V_liq = self.V_tank - V_gas

        ### NOTE: not simulating m_dot_cond!
        m_dot_cond = 0# solve_m_dot_condensed(sat_surf, gas_state, V_gas)


        # Net Mass Transfer of Liquid and Gas CV
        m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj)
        
        #then solve the height of the gas wall
        h_gas_wall = V_gas / (0.25*np.pi*(self.diam_in**2))
        V_gas_wall = 0.25*np.pi*((self.diam_out**2)-(self.diam_in**2))*h_gas_wall
        m_gas_wall = self.rho_wall*V_gas_wall

        h_liq_wall = self.height_tank - h_gas_wall
        V_liq_wall = 0.25*np.pi*((self.diam_out**2)-(self.diam_in**2))*h_liq_wall
        m_liq_wall = self.rho_wall*V_liq_wall

        # NOTE: ADIABATIC FOR DEBUGGING
        # Heat transfer (4) [natural convection] from liq wall to liq
        Q_dot_liq_wall_to_liq = solve_Q_dot_natural_convection_liq(rho_liq, T_wall_liq, T_liq, T_liq, P_tank, 0.021, 0.4, h_liq_wall, (np.pi*self.diam_in*h_liq_wall + self.Tank_Inner_Area), "N2O" ) #relative to liq cv       
        # Heat transfer (5) [natural convection] from gas to gas wall
        Q_dot_gas_wall_to_gas = solve_Q_dot_natural_convection_gas(rho_gas, T_wall_gas, T_gas, T_gas, P_tank, 0.021, 0.4, h_gas_wall, (np.pi*self.diam_in*h_gas_wall + self.Tank_Inner_Area), "N2O" ) #relative to gas cv
       

        Q_dot_liq = Q_dot_liq_wall_to_liq - Q_dot_liq_to_sat_surf
        Q_dot_gas = Q_dot_gas_wall_to_gas + Q_dot_sat_surf_to_gas



        #NOTE: USE ambient properties for air, T_2 will be respective wall temperature (RK var)
        # (6) [natural convection] from atm to liq wall
        Q_dot_atm_to_liq_wall = solve_Q_dot_natural_convection_gas(self.rho_atm, self.T_atm, T_wall_liq, self.T_atm, self.P_atm, 0.59, 0.25, self.height_tank, (np.pi*self.diam_out*h_liq_wall + self.Tank_Inner_Area), "Air") #relative to wall_liq cv
        # (7) [natural convection] from atm to gas wall
        Q_dot_atm_to_gas_wall = solve_Q_dot_natural_convection_gas(self.rho_atm, self.T_atm, T_wall_gas, self.T_atm, self.P_atm, 0.59, 0.25, self.height_tank, (np.pi*self.diam_out*h_gas_wall + self.Tank_Inner_Area), "Air") #relative to wall_gas cv
        # (8) [conduction] from liq wall to gas wall 
        Q_dot_liq_wall_to_gas_wall = solve_Q_dot_conduction( (T_wall_liq-T_wall_gas), self.height_tank, self.k_w, self.Tank_xsec_Area) #relative to wall_liq cv

        # Iteratively solve change in CV Volume
        V_dot_liq = self.V_dot_liq_prev #initial guess for dV_dt_liq
        
        sol = root_scalar(
            lambda V_dot: P_dot_error(
                V_dot, liq_state, sat_surf, gas_state,
                P_tank, m_liq, m_gas, V_liq, V_gas,
                m_dot_inj, m_dot_evap, m_dot_cond,
                Q_dot_liq, Q_dot_gas
            ),
            method="secant",
            x0=V_dot_liq,
            x1=V_dot_liq * 0.99 if V_dot_liq != 0 else 1e-8,  # perturbation
            xtol=self.P_dot_err_tolerance,
            maxiter=100
        )

        if not sol.converged:
            raise RuntimeError("P_dot_error secant solver did not converge")

        V_dot_liq = sol.root        
        
        #while np.abs(P_dot_error(V_dot_liq, liq_state, sat_surf, gas_state, P_tank, m_liq, m_gas, V_liq, V_gas, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas)) > self.P_dot_err_tolerance:
        #    V_dot_liq = secant((lambda V_dot: P_dot_error(V_dot, liq_state, sat_surf, gas_state, P_tank, m_liq, m_gas, V_liq, V_gas, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas)), V_dot_liq)
        


        ### Wall nodes:
        
        height_dot = V_dot_liq / (0.25*np.pi*(self.diam_in**2))

        m_dot_liq_wall = self.rho_wall*(0.25*np.pi*height_dot*((self.diam_out**2)-(self.diam_in**2)))  #BUG: this might be a bit unstable w runge kutta steps?
        m_dot_gas_wall = -m_dot_liq_wall

        #NOTE: for T_dot_wall_liq:  IN: (6)      OUT: (4) AND (8)
        T_dot_wall_liq = ( Q_dot_atm_to_liq_wall - Q_dot_liq_wall_to_liq - Q_dot_liq_wall_to_gas_wall + m_dot_liq_wall*CW*(T_wall_liq - T_wall_gas) ) / (CW*m_liq_wall)

        #NOTE: for T_dot_wall_vap:  IN: (7) and (8)      OUT: (5)
        T_dot_wall_gas = ( Q_dot_atm_to_gas_wall - Q_dot_gas_wall_to_gas + Q_dot_liq_wall_to_gas_wall + m_dot_gas_wall*CW*(T_wall_liq - T_wall_gas) ) / (CW*m_gas_wall)

        T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_liq, liq_state, sat_surf, gas_state, P_tank, m_liq, m_gas, V_liq, V_gas, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, True)

        #print("T_dot_liq_gas: ", T_dot_liq, T_dot_gas)

        #print("rk vars: ", P_tank, T_liq, T_gas, m_liq, m_gas)

        return [T_dot_liq, T_dot_gas, m_dot_liq, m_dot_gas, T_dot_wall_liq, T_dot_wall_gas]








    #this one is easy! 
    def system_of_gas_odes(self, t, y, P_cc):

        _, T_gas, _, m_gas, _, T_wall_gas = y  # Unpack state variables

        rho_gas = m_gas/self.V_tank #know V_tank, m_gas so we can find rho_gas!


        # Mass transfer (1) from injector, then mass balance m_dot_gas = m_dot_inj
        #m_dot_gas = spi_model(self.Cd_1, self.A_inj_1, P_tank, P_cc, rho_gas)

        gas_state = SpanWagnerEOS_SingleState(rho_gas, T_gas) # I get this might be clunky but I was having serious speed issues so I tried a lightweight pressure sol in the iterative method. Yes, there are better ways to setup this entire program in general, my bad, im a lot older and more pythonic than when i started this script

        # Mass transfer (1) from injector
        
        self.state["P_1"] = gas_state.P
        self.state["P_2"] = P_cc
        self.state["rho_1"] = rho_gas
        self.state["h_1"] = gas_state.h
        self.state["T_1"] = T_gas
        
        ## NOTE: how to integrate this? - should it just be 0 ???
        self.state["x_1"] = 0# self.x_tank

        m_dot_gas = -self.injector.m_dot(self.state)

        # Heat transfer (5) [natural convection] from gas to gas wall
        Q_dot_gas_wall_to_gas = 0# solve_Q_dot_natural_convection_gas(rho_gas, T_wall_gas, T_gas, T_gas, gas_state.P, 0.021, 0.4, self.height_tank, (np.pi*self.diam_in*self.height_tank), "N2O" ) #relative to gas cv

        U_dot_gas = m_dot_gas*gas_state.h + Q_dot_gas_wall_to_gas

        d_rho_dt_gas = m_dot_gas/self.V_tank

        T_dot_gas = (1/gas_state.cv)*( (1/m_gas) * (U_dot_gas - (gas_state.u * m_dot_gas)) - (gas_state.du_drho_const_T * d_rho_dt_gas) )

        Q_dot_atm_to_gas_wall = 0 # solve_Q_dot_natural_convection_gas(self.rho_atm, self.T_atm, T_wall_gas, self.T_atm, self.P_atm, 0.59, 0.25, self.height_tank, (np.pi*self.diam_in*self.height_tank), "Air") #relative to wall_gas cv

        T_dot_wall_gas = Q_dot_atm_to_gas_wall/(CW*self.m_tank)

        #print("rk vars: ", gas_state.P, T_gas, m_gas)

        return [0.0, T_dot_gas, 0.0, m_dot_gas, 0.0, T_dot_wall_gas]


    def tank_ode_system(self, t, y, P_cc):
        """
        Wrapper that chooses liquid-phase or vapor-only ODEs.
        """
        _, _, m_liq, _, _, _ = y

        #print("m_liq in non eq tank: ", t, m_liq)

        if m_liq > 0:
            return self.system_of_liq_odes(t, y, P_cc)
        else:
            #we are in the vapor phase (liquid fully drained from tank)
            return self.system_of_gas_odes(t, y, P_cc)












    def inst(self, P_cc):

        t = 0
        y0 = [self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas]
        y_new = rk4_step(self.tank_ode_system, 0.0, y0, self.timestep, P_cc)

        self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas = y_new

        print(f"t = {self.t:.3f} s  |  non eq: {self.T_liq:.3f}, {self.T_gas:.3f}, {self.m_liq:.3f}, {self.m_gas:.3f}, {self.T_wall_liq:.3f}, {self.T_wall_gas:.3f} ")

        # (4) iteratively solve P_tank to update thermodynamic properties in each node
        self.rho_liq, self.rho_gas, self.P_tank = solve_thermo_params(self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.rho_liq_prev, self.rho_gas_prev, self.V_tank)

        #update stored vals for RK est and volumes
        self.rho_liq_prev = self.rho_liq
        self.rho_gas_prev = self.rho_gas
        self.V_dot_liq_prev = (self.V_liq - self.m_liq/self.rho_liq)/self.timestep
        self.V_liq = self.m_liq/self.rho_liq
        self.V_gas = self.V_tank - self.V_liq

        T_sat = self.T_gas #override in liq phase

        if self.m_liq >= 0.0:
            liq_state = SpanWagnerEOS_SingleState(self.rho_liq, self.T_liq)

            self.state["P_1"] = self.P_tank
            self.state["P_2"] = P_cc
            self.state["rho_1"] = self.rho_liq
            self.state["h_1"] = liq_state.h
            self.state["T_1"] = self.T_liq

            sat_surf = SpanWagnerEOS_EquilibriumPhase(None, self.P_tank)
            T_sat = sat_surf.T
        else:
            gas_state = SpanWagnerEOS_SingleState(self.rho_gas, self.T_gas)

            self.state["P_1"] = gas_state.P
            self.state["P_2"] = P_cc
            self.state["rho_1"] = self.rho_gas
            self.state["h_1"] = gas_state.h
            self.state["T_1"] = self.T_gas

        

        m_dot= self.injector.m_dot(self.state) #this one pos

        # Apply smoothing like equilibrium tank
        if self.t <= (self.timestep*2):   # first step
            m_dot = 0.5 * m_dot
        else:
            m_dot = 0.5 * m_dot + 0.5 * self.m_dot_ox_prev
        # Store previous for next call
        self.m_dot_ox_prev = m_dot

        # Advance simulation clock
        self.t += self.timestep

        if self.m_liq > 0.0:
            return {"m_ox_tank": (self.m_gas+self.m_liq), "m_dot_ox": m_dot, "P_ox_tank": self.P_tank, "T_liq_ox_tank": self.T_liq, "T_sat_ox_tank": T_sat, "T_gas_ox_tank": self.T_gas, "state": self.state }
        else: #gas phase
            return {"m_ox_tank": (self.m_gas+self.m_liq), "m_dot_ox": m_dot, "P_ox_tank": self.P_tank, "T_liq_ox_tank": self.T_gas, "T_sat_ox_tank": self.T_gas, "T_gas_ox_tank": self.T_gas, "state": self.state }



















"""
# [7] Test case 1 
t = 0
TIMESTEP = 1e-3

P_atm = 1e5 #Pa
T_atm = 286.5 #K
rho_atm = 1.225 #kg/m^3


### Karabeyoglu test case inputs ~ don't know injector used well!
m_nos = 20 #kg
P_tank = 45e5 #Pa
V_tank = 0.0354 #m^3

# dimensions of test case [7] provided in [8]
diam_in = 0.1905 #m
diam_out = diam_in + (3.18/1000)*2 #m 
rho_wall = 2770 #kg/m^3
k_w = 237 #W/(m K)

Cd_1 = 0.425
A_inj_1 = 0.00003 #m^3 NOTE: GUESS
P_cc = 1.03e6 #Pa

published_time_arr = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5]

published_P_tank_arr = [4.5e6, 4.0e6, 3.75e6, 3.625e6, 3.5e6, 3.375e6, 3.25e6, 3.15e6, 3.05e6, 2.95e6, 2.875e6, 2.625e6, 2.1e6, 1.45e6 ]


#published_m_liq_arr = []
#published_m_gas_arr = []


published_T_liq_arr = [286.5, 286.0, 284.0, 282.5, 280.5, 279.0, 277.5, 275.0, 273.0, 271.0, 269.5, 267.5, None, None]
published_T_gas_arr = [286.5, 282.5, 280.25, 279.0, 277.5, 276.25, 275.0, 274.0, 272.5, 270.25, 269.0, 267.5, None, None]
published_T_sat_arr = [286.5, 282.5, 280.25, 279.0, 277.5, 276.25, 275.0, 274.0, 272.5, 272.0, 271.0, 270.0, None, None] #NOTE: T_sat_arr slightly higher for the first 0.5s but timesteps not high enough fidelity for that

"""

#Tomacz run tank inputs
"""
m_nos = 0.180 #kg
P_tank = 5.2e6 #Pa

diam_in = (40e-3) #m
diam_out = (47.5e-3) #m #bit of a guess
V_tank = 0.25*np.pi*( (40e-3)**2 )*(220e-3)

rho_wall = 2770 #kg/m^3
k_w = 237 #W/(m K)

Cd_1 = 0.45 #SPI discharge coeff for inj #2
A_inj_1 = 0.25*np.pi*((1.5e-3)**2) #m^2
P_cc = P_atm
"""

"""
#def __init__(self, timestep, m_nos, P_tank, P_cc, P_atm, T_atm, rho_atm, V_tank, diam_out, diam_in, rho_wall, k_w, volume_err_tol, P_dot_err_tol, injector):
volume_err_tolerance = 1e-8
P_dot_err_tolerance = 10
#TODO: update: tank = non_equilibrium_tank_model(TIMESTEP, m_nos, P_tank, P_cc, P_atm, T_atm, rho_atm, V_tank, diam_out, diam_in, rho_wall, k_w, Cd_1, A_inj_1, volume_err_tolerance, P_dot_err_tolerance)

time_arr = []
P_tank_arr = []
m_tank_arr = []
m_liq_arr = []
m_gas_arr = []

T_liq_arr = []
T_gas_arr = []
T_sat_arr = []

P_sat_liq_arr = []
P_sat_gas_arr = []
cool_P_sat_liq_arr = []
cool_P_sat_gas_arr = []

T_liq_wall_arr = []
T_gas_wall_arr = []

m_inj_arr = []
U_liq_arr = []
U_gas_arr = []
U_inj_arr = []

V_gas_arr = []
V_liq_arr = []
rho_gas_arr = []
rho_liq_arr = []


try:
    #start_time = time.time()  # Start timer

    while(t < 6000*TIMESTEP): #3000*TIMESTEP
        
        tank.inst(P_cc)
        t+=TIMESTEP 
        LOOKUP_TIME = t

        time_arr.append(t)
        #print("\n next timestep \n")
        
        P_tank_arr.append(tank.P_tank)
        m_liq_arr.append(tank.m_liq)
        m_gas_arr.append(tank.m_gas)




        T_liq_arr.append(tank.T_liq)
        try:
            sat_surf = SpanWagnerEOS_EquilibriumPhase(None, tank.P_tank)
            T_sat_arr.append(sat_surf.T)
        except Exception:
            T_sat_arr.append(0)
        T_gas_arr.append(tank.T_gas)

        T_liq_wall_arr.append(tank.T_wall_liq)
        T_gas_wall_arr.append(tank.T_wall_gas)


        print(f"at t = {t}, final y: ",tank.T_liq, tank.T_gas, tank.m_liq, tank.m_gas, tank.T_wall_liq, tank.T_wall_gas, "\n")

    #end_time = time.time()  # End timer
    #print(f"\nTotal simulation time: {end_time - start_time:.3f} seconds")

except Exception as e:
    traceback.print_exc()




plt.subplot(1,3,1)
plt.scatter(time_arr,P_tank_arr,label = "model_tank", color = "b" )

plt.plot(published_time_arr,published_P_tank_arr,label = "published_tank", color = "dodgerblue" )
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Pa)')
plt.title('Pressure vs. Time')
plt.legend()
plt.grid(True)

plt.subplot(1,3,2)
plt.scatter(time_arr,m_liq_arr, label = "model_liquid", color = "b" )
plt.scatter(time_arr,m_gas_arr, label = "model_gas", color = "darkorange" )

plt.plot(published_time_arr,published_m_liq_arr, label = "published_liquid", color = "dodgerblue" )
plt.plot(published_time_arr,published_m_gas_arr, label = "published_gas", color = "gold" )
plt.xlabel('Time (s)')
plt.ylabel('Mass (kg)')
plt.title('Mass vs. Time')
plt.legend()
plt.grid(True)

plt.subplot(1,3,3)
plt.scatter(time_arr,T_liq_arr, label = "model_liquid", color = "b" )
plt.scatter(time_arr,T_gas_arr, label = "model_gas", color = "darkorange" )
plt.scatter(time_arr,T_sat_arr, label = "model_T_sat", color =  "g" )

plt.plot(published_time_arr,published_T_liq_arr, label = "published_liquid", color = "dodgerblue" )
plt.plot(published_time_arr,published_T_gas_arr, label = "published_gas", color = "gold" )
plt.plot(published_time_arr,published_T_sat_arr, label = "published_T_sat", color = "limegreen" )

plt.scatter(time_arr,T_liq_wall_arr, label = "model_WALL liquid", color = "r" )
plt.scatter(time_arr,T_gas_wall_arr, label = "model_WALL gas", color = "mediumorchid" )
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.title('Temperature vs. Time')
plt.legend()
plt.grid(True)
plt.show()
"""