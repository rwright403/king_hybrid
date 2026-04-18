import numpy as np
from dataclasses import dataclass, field
from scipy.optimize import root, root_scalar, brentq
from src.models.ox_tank._base import BaseTank
from thermo import Chemical
import matplotlib.pyplot as plt
import traceback

import time
from src.models._thermo.convection_heat_transfer import *
from src.models._thermo.conduction_heat_transfer import *

from src.models._thermo.n2o_thermo_span_wagner_class import SpanWagnerEOS_SingleState, SpanWagnerEOS_EquilibriumPhase, lightweight_span_wagner_eos_pressure, lightweight_span_wagner_eos_cp, lightweight_span_wagner_eos_d_rho_dT_P

from src.utils.numerical_methods import rk4_step


from rocketprops.rocket_prop import get_prop
from rocketprops.line_supt import calc_line_vel_dp

#thermo table import
from src.models._thermo.cv_lookup_tables.n2o_cv_gas_lookup import N2OCVGasTable
cv_gas_table = N2OCVGasTable()
from src.models._thermo.cv_lookup_tables.n2o_cv_liq_lookup import N2OCVLiqTable
cv_liq_table = N2OCVLiqTable()
from src.models._thermo.dP_drho_const_T_lookup_tables.n2o_dP_drho_const_T_gas_lookup import N2ODPDRHOTGasTable
dP_drho_const_T_gas_table = N2ODPDRHOTGasTable()
from src.models._thermo.dP_drho_const_T_lookup_tables.n2o_dP_drho_const_T_liq_lookup import N2ODPDRHOTLiqTable
dP_drho_const_T_liq_table = N2ODPDRHOTLiqTable()
from src.models._thermo.dP_dT_const_rho_lookup_tables.n2o_dP_dT_const_rho_gas_lookup import N2ODPDTRHOGasTable
dP_dT_const_rho_gas_table = N2ODPDTRHOGasTable()
from src.models._thermo.dP_dT_const_rho_lookup_tables.n2o_dP_dT_const_rho_liq_lookup import N2ODPDTRHOLiqTable
dP_dT_const_rho_liq_table = N2ODPDTRHOLiqTable()
from src.models._thermo.du_drho_const_T_lookup_tables.n2o_du_drho_const_T_gas_lookup import N2ODUDRHOTGasTable
du_drho_const_T_gas_table = N2ODUDRHOTGasTable()
from src.models._thermo.du_drho_const_T_lookup_tables.n2o_du_drho_const_T_liq_lookup import N2ODUDRHOTLiqTable
du_drho_const_T_liq_table = N2ODUDRHOTLiqTable()
from src.models._thermo.h_lookup_tables.n2o_h_gas_lookup import N2OHGasTable
h_gas_table = N2OHGasTable()
from src.models._thermo.h_lookup_tables.n2o_h_liq_lookup import N2OHLiqTable
h_liq_table = N2OHLiqTable()
from src.models._thermo.P_lookup_tables.n2o_P_gas_lookup import N2OPGasTable
P_gas_table = N2OPGasTable()
from src.models._thermo.P_lookup_tables.n2o_P_liq_lookup import N2OPLiqTable
P_liq_table = N2OPLiqTable()
from src.models._thermo.rho_sat_lookup_tables.n2o_rho_sat_gas_lookup import N2ORhoSatGasTable
rho_sat_gas_table = N2ORhoSatGasTable()
from src.models._thermo.rho_sat_lookup_tables.n2o_rho_sat_liq_lookup import N2ORhoSatLiqTable
rho_sat_liq_table = N2ORhoSatLiqTable()
from src.models._thermo.T_sat_lookup_table.n2o_T_sat_lookup import N2OSatTemperatureTable
T_sat_table = N2OSatTemperatureTable()
from src.models._thermo.P_sat_lookup_table.n2o_P_sat_lookup import N2OSatPressureTable
P_sat_table = N2OSatPressureTable()
from src.models._thermo.u_lookup_tables.n2o_u_gas_lookup import N2OUGasTable
u_gas_table = N2OUGasTable()
from src.models._thermo.u_lookup_tables.n2o_u_liq_lookup import N2OULiqTable
u_liq_table = N2OULiqTable()

# Global Constants:
R_U = 8.31446 #J/(mol K)
T_REF = 298.15 #K

n2o_g = Chemical('nitrous oxide', T=T_REF)
PC = n2o_g.Pc
TC = n2o_g.Tc
OMEGA = n2o_g.omega

MW = 0.0440128 #n2o.MW in g/mol --> converted to kg/mol
print("MW: ", MW)
KAPPA = 0.37464 + 1.5422*n2o_g.omega - 0.26992*n2o_g.omega**2
b = 0.07780*(R_U*TC/PC)
g = 9.81 #m/s^2

CW = 896 #J/(kg K) ~ this is the specific heat capacity of the wall material, for Al 6061 from MATWEB: https://www.matweb.com/search/datasheet.aspx?MatGUID=b8d536e0b9b54bd7b69e4124d8f1d20a&ckck=1

E = 2.1e4 #1.3e3 #2.1e4 #NOTE: EMPIRICAL FACTOR E for high mass flow effects to scale Q_liq_to_sat_surf and m_dot_evap. This was obtained from [7], which claims E is an oxidizer property. This is contradicted by [8] which claims E is a function of the tank as a system



def spi_model(Cd_spi, A, P_1, P_2, rho):

    m_dot_spi = -Cd_spi*A*np.sqrt(2*rho*(P_1-P_2))

    return m_dot_spi


def solve_m_dot_evap(h_sat_gas, h_sat_liq, h_liq, Q_dot_liq_to_sat_surf, Q_dot_sat_surf_to_gas):
    m_dot_evap = 0

    if (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) > 0:

        m_dot_evap = (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) / ( (h_sat_gas-h_sat_liq) + (h_sat_liq - h_liq)  ) 

    return m_dot_evap



def solve_m_dot_condensed(T_gas, P_tank, V_gas):
    m_dot_cond = 0

    #p_sat_gas = CP.PropsSI('P', 'T', gas_state.T, 'Q', 1, "N2O")
    #TODO: replace ^ with 

    p_sat_gas = P_sat_table.lookup(T_gas)
        
    t = 0.1 #relaxation time

    if (P_tank >= p_sat_gas):
        m_dot_cond = ((P_tank-p_sat_gas)*V_gas*MW)/( R_U*T_gas*t )

    # old
    # if (P_tank >= p_sat_gas):
    #     m_dot_cond = ((gas_state.P-p_sat_gas)*V_gas*MW)/( R_U*gas_state.T*t )   #NOTE EDIT DENOM FOR TESTING, OLD FOR REF: ( preos_g.Z_g*(R_U/MW)*T_gas*(TIMESTEP) )


    #NOTE: CONDENSATION FROM PREOS
    #if p_tank >= p_sat_gas, then condensation to enforce equilibrium
    
    #print(f"m_dot_cond: {m_dot_cond:.5f}")

    return m_dot_cond



def solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj, m_dot_vent, m_dot_liq_in):
    m_dot_liq = -m_dot_evap + m_dot_cond + m_dot_liq_in + m_dot_inj #note inj signed (-) so added
    m_dot_gas =  (m_dot_vent) + m_dot_evap - m_dot_cond #only m_dot_vent sign convention rel to gas node: for remaining terms, convert sign convention from liq cv to gas cv

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

def P_error(m_gas, V_tank_remaining, T_atm, P_atm):

    rho_gas = m_gas/V_tank_remaining
    P_gas = lightweight_span_wagner_eos_pressure(rho_gas, T_atm)

    P_residual = P_atm - P_gas 

    return P_residual





def single_solve_T_dot_liq_gas(V_dot_liq, P_tank, m_liq, m_gas, V_liq, V_gas, cv_liq, cv_gas, h_gas, h_sat_gas, h_sat_liq, h_liq, u_liq, u_gas, du_drho_const_T_liq, du_drho_const_T_gas, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, m_dot_vent, m_dot_liq_in, h_liq_in, debug_mode):

    m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj, m_dot_vent, m_dot_liq_in)

    V_dot_gas = -V_dot_liq

    d_rho_dt_liq = (1/V_liq)*m_dot_liq -(m_liq/(V_liq**2))*V_dot_liq
    d_rho_dt_gas = (1/V_gas)*m_dot_gas -(m_gas/(V_gas**2))*V_dot_gas


    Q_dot_cond_release = m_dot_cond * (h_gas - h_sat_liq)

    U_dot_liq = (
        m_dot_liq_in*h_liq_in
        + m_dot_inj*h_liq
        - m_dot_evap*h_liq
        + m_dot_cond*h_sat_liq
        - P_tank*V_dot_liq
        + Q_dot_liq
    )

    U_dot_gas = (
        m_dot_vent*h_gas         
        + m_dot_evap*h_sat_gas
        - m_dot_cond*h_sat_liq
        - P_tank*V_dot_gas
        + Q_dot_gas
        + Q_dot_cond_release
    )


    T_dot_liq = (1/cv_liq)*( (1/m_liq) * (U_dot_liq - (u_liq * m_dot_liq)) - (du_drho_const_T_liq * d_rho_dt_liq) )
    T_dot_gas = (1/cv_gas)*( (1/m_gas) * (U_dot_gas - (u_gas * m_dot_gas)) - (du_drho_const_T_gas * d_rho_dt_gas) )


    if debug_mode == True:
        a=1
        #print("T_dot_liq, T_dot_gas: ", T_dot_liq, T_dot_gas)
        #print("T_dot_gas: ", T_dot_gas, U_dot_gas,  m_dot_evap*(sat_surf.h_sat_gas), - m_dot_cond*(gas_state.h), - P_tank*V_dot_gas, + Q_dot_gas)
        #print("T_dot_gas: ", T_dot_gas, (1/gas_state.cv)*(1/m_gas)*(U_dot_gas - (gas_state.u * m_dot_gas)) , -(1/gas_state.cv)*(gas_state.du_drho_const_T * d_rho_dt_gas) )


    return T_dot_liq, T_dot_gas 

def P_dot_error(V_dot_guess, P_tank, m_liq, m_gas, V_liq, V_gas, cv_liq, cv_gas, h_gas, h_sat_gas, h_sat_liq, h_liq, u_liq, u_gas, du_drho_const_T_liq, du_drho_const_T_gas, dP_dT_const_rho_liq, dP_dT_const_rho_gas, dP_drho_const_T_liq, dP_drho_const_T_gas, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, m_dot_vent, m_dot_liq_in, h_liq_in):    

    V_dot_gas = -V_dot_guess #guessing for liquid

    m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj, m_dot_vent, m_dot_liq_in) #or put this outside and pass in?


    d_rho_dt_liq = (1/V_liq)*m_dot_liq - (m_liq/(V_liq**2))*V_dot_guess
    d_rho_dt_gas = (1/V_gas)*m_dot_gas - (m_gas/(V_gas**2))*V_dot_gas

    T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_guess, P_tank, m_liq, m_gas, V_liq, V_gas, cv_liq, cv_gas, h_gas, h_sat_gas, h_sat_liq, h_liq, u_liq, u_gas, du_drho_const_T_liq, du_drho_const_T_gas, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, m_dot_vent, m_dot_liq_in, h_liq_in, False)

    P_dot_liq = dP_dT_const_rho_liq*T_dot_liq + dP_drho_const_T_liq*d_rho_dt_liq

    P_dot_gas = dP_dT_const_rho_gas*T_dot_gas + dP_drho_const_T_gas*d_rho_dt_gas

    return P_dot_liq - P_dot_gas


class non_equilibrium_tank_model():
    def __init__(self, timestep, name, m_ox, P_tank, P_atm, T_atm, rho_atm, V_tank, diam_out, diam_in, rho_wall, k_w, Cd_inj, A_inj, Cd_vent, A_vent, volume_err_tol, P_dot_err_tol):
        self.name = name
        
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

        if self.name == "bottle":
            # enough nitrous such that there is both liq and gas in tank

            #print(f"{self.name} init cond printing temps: | T_eq, T_amb, diff: {initial_eq_state.T:.3f}, {T_atm:.3f}, {(initial_eq_state.T-T_atm):.3f}")

            self.rho_liq = initial_eq_state.rho_sat_liq
            self.rho_gas = initial_eq_state.rho_sat_gas
            
            rho_bulk_tank = self.m_ox/V_tank
            x_tank = ( (1/rho_bulk_tank)-(1/self.rho_liq) ) / ( (1/self.rho_gas)-(1/self.rho_liq) )


            #gas cv setup
            self.T_gas = initial_eq_state.T # bad --> T_atm
            self.m_gas = x_tank*self.m_ox
            v_gas = 1/self.rho_gas
            self.V_gas = v_gas*self.m_gas

            #liquid cv
            self.T_liq = initial_eq_state.T # bad --> T_atm
            self.m_liq = self.m_ox-self.m_gas
            v_liq = 1/self.rho_liq
            self.V_liq = v_liq*self.m_liq

            #print(f"bottle m_liq: {self.m_liq}, m_gas: {self.m_gas:.3f}, x_tank: {x_tank:.3f}, rho_bulk_tank: {rho_bulk_tank:.3f}")


        elif self.name == "run_tank":
            # instead of using air and mixtures assuming tank is filled w nitrous at P_atm, T_atm for run tank init cond. and a small bit of seed mass

            #NOTE: liq and gas NOT at thermal equilibrium



            #NOTE: liq node at quality 0 sat phase
            self.rho_liq = initial_eq_state.rho_sat_liq

            self.T_liq = initial_eq_state.T 
            self.m_liq = 1e-2 #kg
            v_liq = 1/self.rho_liq
            self.V_liq = v_liq*self.m_liq
            

            #NOTE: gas solved at P_atm, T_atm, don't know mass
            self.T_gas = T_atm


            m_gas_guess = 0.5 #kg NOTE: initial guess for secant
            V_tank_remaining = self.V_tank - self.V_liq

            sol = root_scalar(
                lambda m_gas_guess: P_error(m_gas_guess, V_tank_remaining, self.T_atm, self.P_atm),
                method="secant",
                x0=m_gas_guess,
                x1=m_gas_guess * 0.99,
                xtol=1e-9,
                maxiter=100
            )

            if not sol.converged:
                raise RuntimeError("P_dot_error secant solver did not converge")

            self.m_gas = sol.root


            self.V_gas = V_tank_remaining
            self.rho_gas = self.m_gas / self.V_gas

            #print("m_gas, rho_gas: ", self.m_gas, self.rho_gas)

            #print(f"run_tank: T_liq = {self.T_liq:.3f} T_gas = {self.T_gas:.3f} m_liq = {self.m_liq} m_gas = {self.m_gas:.3f}, rho_liq = {self.rho_liq:.3f}, rho_gas = {self.rho_gas:.3f} ")

        else:
            raise NameError("wrong name - doesn't correspond to an init cond, options are: 'run_tank' or 'bottle'\n")


        self.V_tank = self.V_liq+self.V_gas # resolving just to get rid of unstability "what are you going to do if the aluminum is too small? water it? give it sunlight? let it grow?"
        self.height_tank = self.V_tank/(0.25*np.pi*(diam_in**2))

        # NOTE: ASSUME TANK WALL AT T_ATM, MAKE SURE THIS MATCHES REAL WORLD INITIAL CONDITIONS
        self.T_wall_gas = T_atm 
        self.T_wall_liq = T_atm


        self.rho_liq_prev = self.rho_liq
        self.rho_gas_prev = self.rho_gas

        self.V_dot_liq_prev = 1E-8

        self.t = 0.0
        # self.m_dot_vent_prev
        #TODO: CHANGE BELOW TO INJ
        self.m_dot_inj_out_prev = 0.0


        ### mass flow consts:
        self.Cd_inj = Cd_inj
        self.A_inj = A_inj
        self.Cd_vent = Cd_vent
        self.A_vent = A_vent


        self.m_tank = 0.25*np.pi*(self.diam_out**2-self.diam_in**2) * self.height_tank

        # x tank not in both scope! print("init non eq tank: ", x_tank, self.V_tank, self.rho_liq, rho_bulk_tank, self.rho_gas)
   


    #TODO: create parms vector
    def rhs(self, t, y, P_cc, m_dot_liq_in, h_liq_in):
        

        #print(f"\nname: {self.name}")

        T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas = y  # Unpack state variables

        #print(f"in {self.name}, T_liq = {T_liq:.3f}, T_gas = {T_gas:.3f}, m_liq = {m_liq}, m_gas = {m_gas:.3f} ")

        ### Solve thermo parameters! - old according to [8]
        
        rho_liq, rho_gas, P_tank = solve_thermo_params(T_liq, T_gas, m_liq, m_gas, self.rho_liq_prev, self.rho_gas_prev, self.V_tank)


        #print("after root sol: ", T_liq, T_gas, rho_liq, rho_gas, P_tank)




        cv_gas = cv_gas_table.lookup(rho_gas, T_gas)
        dP_drho_const_T_gas = dP_drho_const_T_gas_table.lookup(rho_gas, T_gas)
        dP_dT_const_rho_gas = dP_dT_const_rho_gas_table.lookup(rho_gas, T_gas)
        du_drho_const_T_gas = du_drho_const_T_gas_table.lookup(rho_gas, T_gas)
        h_gas = h_gas_table.lookup(rho_gas, T_gas)
        u_gas = u_gas_table.lookup(rho_gas, T_gas)
        
        
        
        
        T_sat = T_sat_table.lookup(P_tank)
        rho_sat_gas = rho_sat_gas_table.lookup(P_tank)
        rho_sat_liq = rho_sat_liq_table.lookup(P_tank)
        h_sat_gas = h_gas_table.lookup(rho_sat_gas, T_sat)
        h_sat_liq = h_liq_table.lookup(rho_sat_liq, T_sat)





        cv_liq = cv_liq_table.lookup(rho_liq, T_liq)
        dP_drho_const_T_liq = dP_drho_const_T_liq_table.lookup(rho_liq, T_liq)
        dP_dT_const_rho_liq = dP_dT_const_rho_liq_table.lookup(rho_liq, T_liq)
        du_drho_const_T_liq = du_drho_const_T_liq_table.lookup(rho_liq, T_liq)
        h_liq = h_liq_table.lookup(rho_liq, T_liq)
        u_liq = u_liq_table.lookup(rho_liq, T_liq)



        # Mass transfer (1) from injector
        m_dot_inj = spi_model(self.Cd_inj, self.A_inj, P_tank, P_cc, rho_liq)

        # Also mass transfer from vent
        m_dot_vent = spi_model(self.Cd_vent, self.A_vent, P_tank, P_cc, rho_gas)

        

        # Apply mass flow rate smoothing like equilibrium tank
        if self.t <= (self.timestep*2):   # first step
            m_dot_inj = 0.5 * m_dot_inj
            #TODO? m_dot_vent = 0#  0.5 * m_dot_vent
        else:
            m_dot_inj = 0.5 * m_dot_inj + 0.5 * (-self.m_dot_inj_out_prev)
            #TODO? m_dot_vent = 0.5* m_dot_vent + 0.5 * (-self.m_dot_vent_prev)




        
        # Heat transfer (2) from saturated surface to gas                       (T_1, T_2, P_tank, rho_2, c, n, tank_diam, fluid)
        # L = tank inner diam , Area of circle x section
        # OLD: T_film_gas = ((sat_surf.T + T_gas)/2 )
        # OLD: Q_dot_sat_surf_to_gas = solve_Q_dot_natural_convection_gas(rho_gas, sat_surf.T, T_gas, T_film_gas, P_tank, 0.15, 0.333, self.diam_in, self.Tank_Inner_Area, "N2O" ) #relative to gas cv
        T_film_gas = ((T_sat + T_gas)/2 )
        Q_dot_sat_surf_to_gas = solve_Q_dot_natural_convection_gas(rho_gas, T_sat, T_gas, T_film_gas, P_tank, 0.15, 0.333, self.diam_in, self.Tank_Inner_Area, "N2O" ) #relative to gas cv
        
        # Heat transfer (3)  from liq to saturated surface (sat surface assumed to be a liquid with quality 0)
        # OLD: T_film_liq = ((sat_surf.T + T_liq)/2 )
        # OLD: Q_dot_liq_to_sat_surf = (E)*solve_Q_dot_natural_convection_liq(rho_liq, T_liq, sat_surf.T, T_film_liq, P_tank, 0.15, 0.333, self.diam_in, self.Tank_Inner_Area, "N2O" ) #relative to liq cv
        T_film_liq = ((T_sat + T_liq)/2 )
        Q_dot_liq_to_sat_surf = (E)*solve_Q_dot_natural_convection_liq(rho_liq, T_liq, T_sat, T_film_liq, P_tank, 0.15, 0.333, self.diam_in, self.Tank_Inner_Area, "N2O" ) #relative to liq cv
        #NOTE:CORRECTION FACTOR for nitrous oxide heat transfer E = (E) to account for blowing as per [7],[8]

        # Mass transfer (3) by evaporation 
        m_dot_evap = solve_m_dot_evap(h_sat_gas, h_sat_liq, h_liq, Q_dot_liq_to_sat_surf, Q_dot_sat_surf_to_gas)

        # Mass transfer (2) by condensation
        V_gas = m_gas/rho_gas
        V_liq = self.V_tank - V_gas

        m_dot_cond = solve_m_dot_condensed(T_gas, P_tank, V_gas)


        # Net Mass Transfer of Liquid and Gas CV
        m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj, m_dot_vent, m_dot_liq_in)
        
        #then solve the height of the gas wall
        h_gas_wall = V_gas / (0.25*np.pi*(self.diam_in**2))
        V_gas_wall = 0.25*np.pi*((self.diam_out**2)-(self.diam_in**2))*h_gas_wall
        m_gas_wall = self.rho_wall*V_gas_wall

        h_liq_wall = self.height_tank - h_gas_wall
        V_liq_wall = 0.25*np.pi*((self.diam_out**2)-(self.diam_in**2))*h_liq_wall
        m_liq_wall = self.rho_wall*V_liq_wall


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
                V_dot, P_tank, m_liq, m_gas, V_liq, V_gas,
                cv_liq, cv_gas, h_gas, h_sat_gas, h_sat_liq, h_liq, u_liq, u_gas, du_drho_const_T_liq, du_drho_const_T_gas, dP_dT_const_rho_liq, dP_dT_const_rho_gas, dP_drho_const_T_liq, dP_drho_const_T_gas,
                m_dot_inj, m_dot_evap, m_dot_cond,
                Q_dot_liq, Q_dot_gas, m_dot_vent, m_dot_liq_in, h_liq_in
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



        ### Wall nodes:
        
        height_dot = V_dot_liq / (0.25*np.pi*(self.diam_in**2))

        m_dot_liq_wall = self.rho_wall*(0.25*np.pi*height_dot*((self.diam_out**2)-(self.diam_in**2)))  #BUG: this might be a bit unstable w runge kutta steps?
        m_dot_gas_wall = -m_dot_liq_wall

        #NOTE: for T_dot_wall_liq:  IN: (6)      OUT: (4) AND (8)
        T_dot_wall_liq = ( Q_dot_atm_to_liq_wall - Q_dot_liq_wall_to_liq - Q_dot_liq_wall_to_gas_wall + m_dot_liq_wall*CW*(T_wall_liq - T_wall_gas) ) / (CW*m_liq_wall)

        #NOTE: for T_dot_wall_vap:  IN: (7) and (8)      OUT: (5)
        T_dot_wall_gas = ( Q_dot_atm_to_gas_wall - Q_dot_gas_wall_to_gas + Q_dot_liq_wall_to_gas_wall + m_dot_gas_wall*CW*(T_wall_liq - T_wall_gas) ) / (CW*m_gas_wall)

        T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_liq, P_tank, m_liq, m_gas, V_liq, V_gas, cv_liq, cv_gas, h_gas, h_sat_gas, h_sat_liq, h_liq, u_liq, u_gas, du_drho_const_T_liq, du_drho_const_T_gas, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, m_dot_vent, m_dot_liq_in, h_liq_in, True)




        return [T_dot_liq, T_dot_gas, m_dot_liq, m_dot_gas, T_dot_wall_liq, T_dot_wall_gas]











    #NOTE: this is solver, should do something better!
    def inst(self, P_cc, m_dot_liq_in, h_liq_in):

        t = 0
        y0 = [self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas]
        y_new = rk4_step(self.rhs, 0.0, y0, self.timestep, P_cc, m_dot_liq_in, h_liq_in)

        self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas = y_new

        #print(f"t = {self.t:.3f} s  |  non eq: {self.T_liq:.3f}, {self.T_gas:.3f}, {self.m_liq}, {self.m_gas:.3f}, {self.T_wall_liq:.3f}, {self.T_wall_gas:.3f} ")

        # (4) iteratively solve P_tank to update thermodynamic properties in each node
        self.rho_liq, self.rho_gas, self.P_tank = solve_thermo_params(self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.rho_liq_prev, self.rho_gas_prev, self.V_tank)

        #update stored vals for RK est and volumes
        self.rho_liq_prev = self.rho_liq
        self.rho_gas_prev = self.rho_gas
        self.V_dot_liq_prev = (self.m_liq/self.rho_liq - self.V_liq)/self.timestep
        self.V_liq = self.m_liq/self.rho_liq
        self.V_gas = self.V_tank - self.V_liq



        #old liq_state = SpanWagnerEOS_SingleState(self.rho_liq, self.T_liq)
        #old h_liq = liq_state.h

        h_liq = h_liq_table.lookup(self.rho_liq, self.T_liq)

        #old sat_surf = SpanWagnerEOS_EquilibriumPhase(None, self.P_tank)
        T_sat = T_sat_table.lookup(self.P_tank)



        
        #NOTE: signed (+)  - both inj and vent?
        m_dot_inj_out= -spi_model(self.Cd_inj, self.A_inj, self.P_tank, P_cc, self.rho_liq)


        # Apply smoothing like equilibrium tank
        if self.t <= (self.timestep*2):   # first step
            m_dot_inj_out = 0.5 * m_dot_inj_out
        else:
            m_dot_inj_out = 0.5 * m_dot_inj_out + 0.5 * self.m_dot_inj_out_prev
        # Store previous for next call
        self.m_dot_inj_out_prev = m_dot_inj_out


        # Advance simulation clock
        self.t += self.timestep

        #print("m_dot_inj_out: ", m_dot_inj_out)


        return {"m_ox_tank": (self.m_gas+self.m_liq), "m_dot_inj_out": m_dot_inj_out, "P_ox_tank": self.P_tank, "T_liq_ox_tank": self.T_liq, "T_sat_ox_tank": T_sat, "T_gas_ox_tank": self.T_gas, "h_liq": h_liq}






### Global Consts ###
t = 0
TIMESTEP = 1e-3

P_atm = 1e5 #Pa
T_atm = 286.5 #K
rho_atm = 1.225 #kg/m^3



volume_err_tolerance = 1e-8
P_dot_err_tolerance = 10




### NITROUS BOTTLE ###

botl_P_tank = 45e5 #Pa #TODO: check this across envelope, worst case fill time at T_min

botl_m_nos = 29.4835 #kg #NOTE: 65LB converted to kg
botl_V_tank = 0.1 #m^3 #NOTE: PLACEHOLDER VAL, SHOULD BE: this is for a Q type, can either buy K type or Q type from Linde

botl_diam_in = (0.0254*7) #m #from Linde website
botl_diam_out = botl_diam_in + 2*(0.0254*0.25) #m #guessing, as stated before wall node

# i believe bottle is al but either way i dont believe wall nodes at our scale matter much
botl_rho_wall = 2770 #kg/m^3
botl_k_w = 237 #W/(m K)

botl_Cd_inj = 0.1
botl_A_inj = 0.25*np.pi*(0.0254*(.4))**2

#no draining, placeholder
botl_Cd_vent = 0.65 #using this from RPE 9th ed p279 
botl_A_vent = 0.25*np.pi*(0.0254*(.010))**2 #NOTE: assuming 2'" diam outlet

botl = non_equilibrium_tank_model(TIMESTEP, "bottle", botl_m_nos, botl_P_tank, P_atm, T_atm, rho_atm, botl_V_tank, botl_diam_out, botl_diam_in, botl_rho_wall, botl_k_w, botl_Cd_inj, botl_A_inj, botl_Cd_vent, botl_A_vent, volume_err_tolerance, P_dot_err_tolerance)



### S1 RUN TANK ### - GND % fill
full_pcnt_fill = 0.05
#TODO: UPDATE

s1_m_nos = 0 #kg NOTE: placeholder, solved to obtain density for thermo state at P_atm, T_atm 
s1_P_tank = P_atm #Pa #NOTE: starts as empty, assuming initial condition of tank is filled w nitrous @ P_atm
s1_V_tank = 0.0208 #m^3

s1_diam_in = 0.1282 #m 
s1_diam_out = 0.1413 #m 
s1_rho_wall = 2770 #kg/m^3
s1_k_w = 237 #W/(m K)

s1_Cd_vent = 0.65 #using this from RPE 9th ed p279 
s1_A_vent = 0.25*np.pi*(0.0254*(.002))**2 #NOTE: assuming 

#no draining, placeholder
s1_Cd_inj = 0 #0.57
s1_A_inj = 0 #0.25 * 3.14159 * (42) * (1.5e-3)**2 

"""
### S1 RUN TANK ### - FLIGHT % fill
full_pcnt_fill = 0.9

s1_m_nos = 13.25 #kg
s1_P_tank = P_atm #Pa #NOTE: starts as empty, assuming initial condition of tank is filled w nitrous @ P_atm
s1_V_tank = 0.0208 #m^3

s1_diam_in = 0.1282 #m 
s1_diam_out = 0.1413 #m 
s1_rho_wall = 2770 #kg/m^3
s1_k_w = 237 #W/(m K)

s1_Cd_vent = 0.65 #using this from RPE 9th ed p279 
s1_A_vent = 0.25*np.pi*(0.0254*(.002))**2 #NOTE: assuming 

#no draining, placeholder
s1_Cd_inj = 0.57
s1_A_inj = 0.25 * 3.14159 * (42) * (1.5e-3)**2 
"""


run_tank = non_equilibrium_tank_model(TIMESTEP, "run_tank", s1_m_nos, s1_P_tank, P_atm, T_atm, rho_atm, s1_V_tank, s1_diam_out, s1_diam_in, s1_rho_wall, s1_k_w, s1_Cd_inj, s1_A_inj, s1_Cd_vent, s1_A_vent, volume_err_tolerance, P_dot_err_tolerance)




### Plotting:

@dataclass
class plt_data:
    time_arr: list = field(default_factory=list)
    P_tank_arr: list = field(default_factory=list)
    m_liq_arr: list = field(default_factory=list)
    m_gas_arr: list = field(default_factory=list)
    m_tank_arr: list = field(default_factory=list)

    T_liq_arr: list = field(default_factory=list)
    T_gas_arr: list = field(default_factory=list)
    T_sat_arr: list = field(default_factory=list)

    T_liq_wall_arr: list = field(default_factory=list)
    T_gas_wall_arr: list = field(default_factory=list)



botl_dat = plt_data()
run_tank_dat = plt_data()

def append_dat(t, dat, tank):
    dat.time_arr.append(t)
    dat.P_tank_arr.append(tank.P_tank)
    dat.m_liq_arr.append(tank.m_liq)
    dat.m_gas_arr.append(tank.m_gas)
    #dat.m_tank_arr.append(tank.m_liq+tank.m_gas)

    dat.T_liq_arr.append(tank.T_liq)
    dat.T_gas_arr.append(tank.T_gas)
    try:
        #sat_surf = SpanWagnerEOS_EquilibriumPhase(None, tank.P_tank)
        T_sat = T_sat_table.lookup(tank.P_tank)
        dat.T_sat_arr.append(T_sat)
    except Exception:
        dat.T_sat_arr.append(0)

    dat.T_liq_wall_arr.append(tank.T_wall_liq)
    dat.T_gas_wall_arr.append(tank.T_wall_gas)


def plot_dat(dat, fig_name):
    plt.figure()

    plt.suptitle(fig_name, fontsize=16)

    plt.subplot(1,3,1)
    plt.scatter(dat.time_arr,dat.P_tank_arr,label = "model_tank", color = "b" )

    #plt.plot(published_time_arr,published_P_tank_arr,label = "published_tank", color = "dodgerblue" )
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (Pa)')
    plt.title('Pressure vs. Time')
    plt.legend()
    plt.grid(True)

    plt.subplot(1,3,2)
    plt.scatter(dat.time_arr,dat.m_liq_arr, label = "model_liquid", color = "b" )
    plt.scatter(dat.time_arr,dat.m_gas_arr, label = "model_gas", color = "darkorange" )


    plt.xlabel('Time (s)')
    plt.ylabel('Mass (kg)')
    plt.title('Mass vs. Time')
    plt.legend()
    plt.grid(True)

    plt.subplot(1,3,3)
    plt.scatter(dat.time_arr,dat.T_liq_arr, label = "model_liquid", color = "b" )
    plt.scatter(dat.time_arr,dat.T_gas_arr, label = "model_gas", color = "darkorange" )
    plt.scatter(dat.time_arr,dat.T_sat_arr, label = "model_T_sat", color =  "g" )

    #plt.plot(published_time_arr,published_T_liq_arr, label = "published_liquid", color = "dodgerblue" )
    #plt.plot(published_time_arr,published_T_gas_arr, label = "published_gas", color = "gold" )
    #plt.plot(published_time_arr,published_T_sat_arr, label = "published_T_sat", color = "limegreen" )

    plt.scatter(dat.time_arr,dat.T_liq_wall_arr, label = "model_WALL liquid", color = "r" )
    plt.scatter(dat.time_arr,dat.T_gas_wall_arr, label = "model_WALL gas", color = "mediumorchid" )
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature vs. Time')
    plt.legend()
    plt.grid(True)




try:
    start_time = time.time()  # Start timer

    pcnt_fill = 0.0 #dont need to update
    P_run_tank = run_tank.P_tank

    while( t <=0*TIMESTEP): #pcnt_fill < full_pcnt_fill):
        

        botl_out = botl.inst(P_run_tank, 0, 0)

        #add rocketprops check here

        run_tank.inst(P_atm, botl_out.get("m_dot_inj_out"), botl_out.get("h_liq"))

        #recalc pcnt_fill and P_run_tank for next timestep
        pcnt_fill = run_tank.V_liq/run_tank.V_tank

        P_run_tank = run_tank.P_tank

        #march timestep
        t+=TIMESTEP 
        #LOOKUP_TIME = t

        # logging
        append_dat(t, botl_dat, botl)
        append_dat(t, run_tank_dat, run_tank)

        

        #print(f"in bottle at t = {t}, final y: ", botl.T_liq, botl.T_gas, botl.m_liq, botl.m_gas, botl.T_wall_liq, botl.T_wall_gas, "\n")
        #print(f"in run tank at t = {t}, final y: ", run_tank.T_liq, run_tank.T_gas, run_tank.m_liq, run_tank.m_gas, run_tank.T_wall_liq, run_tank.T_wall_gas, "\n")

    end_time = time.time()  # End timer
    print(f"\nTotal simulation time for one step: {end_time - start_time:.3f} seconds")

except Exception as e:
    traceback.print_exc()

plot_dat(botl_dat, "BOTTLE")
plot_dat(run_tank_dat, "RUN TANK")
plt.show()
