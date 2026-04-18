import numpy as np
from dataclasses import dataclass, field
from scipy.optimize import root, root_scalar, brentq
from src.models.ox_tank._base import BaseTank
import matplotlib.pyplot as plt
import traceback

import time
import cProfile
import pstats

from numba import njit


from src.models._thermo.convection_heat_transfer import *
from src.models._thermo.conduction_heat_transfer import *


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

MW = 0.0440128 #n2o.MW in g/mol --> converted to kg/mol

g = 9.81 #m/s^2

CW = 896 #J/(kg K) ~ this is the specific heat capacity of the wall material, for Al 6061 from MATWEB: https://www.matweb.com/search/datasheet.aspx?MatGUID=b8d536e0b9b54bd7b69e4124d8f1d20a&ckck=1

E = 2.1e4 #1.3e3 #2.1e4 #NOTE: EMPIRICAL FACTOR E for high mass flow effects to scale Q_liq_to_sat_surf and m_dot_evap. This was obtained from [7], which claims E is an oxidizer property. This is contradicted by [8] which claims E is a function of the tank as a system

TIMESTEP = 3e-3

P_atm = 1e5 #Pa
T_atm = 286.5 #K
rho_atm = 1.225 #kg/m^3




def rk4_step(rhs_func, t, y, dt, u, p, guess):
    k1 = rhs_func(t, y, dt, u, p, guess)
    k2 = rhs_func(t + 0.5*dt, y + 0.5*dt*k1, dt, u, p, guess)
    k3 = rhs_func(t + 0.5*dt, y + 0.5*dt*k2, dt, u, p, guess)
    k4 = rhs_func(t + dt, y + dt*k3, dt, u, p, guess)
    return y + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)






def spi_model(Cd_spi, A, P_1, P_2, rho):

    m_dot_spi = -Cd_spi*A*np.sqrt(2*rho*(P_1-P_2))

    return m_dot_spi




#NOTE: using quality at P_tank to distribute m_dot_upstream_in into liq and gas node.
def solve_m_dot_upstream_quality(m_dot_in_upstream, h_in_upstream, P_tank):
    
    T_sat = T_sat_table.lookup(P_tank)
    rho_sat_gas = rho_sat_gas_table.lookup(P_tank)
    rho_sat_liq = rho_sat_liq_table.lookup(P_tank)
    h_sat_gas = h_gas_table.lookup(rho_sat_gas, T_sat)
    h_sat_liq = h_liq_table.lookup(rho_sat_liq, T_sat)
    

    x = (h_in_upstream - h_sat_liq)/(h_sat_gas - h_sat_liq)

    if x >= h_sat_gas:
        x = 1               #x=1 means all the mass in is a vapor
    elif x <= h_sat_liq:
        x = 0                #x=0 means all the mass in is a liquid
    
    m_dot_in_liq = (1-x)*m_dot_in_upstream
    m_dot_in_gas = x*m_dot_in_upstream

    return m_dot_in_liq, m_dot_in_gas


#NOTE: this one enforces the boundary condition between the liq and gas nodes
def solve_m_dot_evap(h_sat_gas, h_sat_liq, h_liq, Q_dot_liq_to_sat_surf, Q_dot_sat_surf_to_gas):
    m_dot_evap = 0


    print("Q_dot sat surfs: ", Q_dot_liq_to_sat_surf, Q_dot_sat_surf_to_gas, -Q_dot_sat_surf_to_gas)

    if (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) > 0 and Q_dot_sat_surf_to_gas>0: #NOTE: should be fixed now

        m_dot_evap = (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) / ( (h_sat_gas-h_sat_liq) + (h_sat_liq - h_liq)  ) 

    return m_dot_evap








def solve_m_dot_condensed(T_gas, P_tank, V_gas):
    m_dot_cond = 0

    p_sat_gas = P_sat_table.lookup(T_gas)
        
    t = 0.1 #relaxation time

    if (P_tank >= p_sat_gas):
        m_dot_cond = ((P_tank-p_sat_gas)*V_gas*MW)/( R_U*T_gas*t )


    return m_dot_cond



def solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj, m_dot_vent, m_dot_liq_in):
    m_dot_liq = -m_dot_evap + m_dot_cond + m_dot_liq_in + m_dot_inj #note inj signed (-) so added
    m_dot_gas =  (m_dot_vent) + m_dot_evap - m_dot_cond #only m_dot_vent sign convention rel to gas node, also signed (-) from inj: for remaining terms, convert sign convention from liq cv to gas cv

    return m_dot_liq, m_dot_gas


def thermo_residuals(rhos, T_liq, T_gas, m_liq, m_gas, V_tank):
    rho_liq, rho_gas = rhos

    P_liq = P_liq_table.lookup(rho_liq, T_liq)
    P_gas = P_gas_table.lookup(rho_gas, T_gas)

    V_est = (m_liq / rho_liq) + (m_gas / rho_gas)


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
    P_tank = P_gas_table.lookup(rho_gas, T_gas)

    V_liq = m_liq/rho_liq
    V_gas = V_tank - V_liq

    return rho_liq, rho_gas, P_tank, V_liq, V_gas



def P_error(m_gas, V_tank_remaining, T_atm, P_atm):

    rho_gas = m_gas/V_tank_remaining
    #P_gas = lightweight_span_wagner_eos_pressure(rho_gas, T_atm)
    P_gas = P_gas_table.lookup(rho_gas, T_atm)

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
        #print("T_dot_gas: ", T_dot_gas, (1/cv_gas)*(1/m_gas)*(U_dot_gas - (u_gas * m_dot_gas)) , -(1/cv_gas)*(du_drho_const_T_gas * d_rho_dt_gas), m_dot_gas )
        
        print(f"U_dot_gas: vent={m_dot_vent*h_gas}, evap={m_dot_evap*h_sat_gas}, cond={-m_dot_cond*h_sat_liq}, PdV={-P_tank*V_dot_gas}, Q_gas={Q_dot_gas}, Q_cond_rel={Q_dot_cond_release}, total={U_dot_gas}")


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


#NEW INIT STATE



def make_tank_geometry_and_consts(diam_in, diam_out, V_tank, rho_wall):

    if diam_in >= diam_out:
        raise ValueError("Tank inner diameter >= Tank outer diameter, this will give nonphysical results so check your inputs!")
    
    Tank_Inner_Area = 0.25*np.pi*(diam_in**2)
    Tank_xsec_Area = 0.25*np.pi*((diam_out**2)-(diam_in**2))

    height_tank = V_tank/Tank_Inner_Area
    m_tank = rho_wall * Tank_xsec_Area * height_tank

    return Tank_Inner_Area, Tank_xsec_Area, height_tank, m_tank






    
def build_bottle_initial_state(P_tank, V_tank, m_ox, T_atm):

    T_sat = T_sat_table.lookup(P_tank)
    rho_gas = rho_sat_gas_table.lookup(P_tank)
    rho_liq = rho_sat_liq_table.lookup(P_tank)

    rho_bulk_tank = m_ox/V_tank
    x_tank = ( (1/rho_bulk_tank)-(1/rho_liq) ) / ( (1/rho_gas)-(1/rho_liq) )


    #gas node setup
    T_gas = T_sat
    m_gas = x_tank*m_ox
    v_gas = 1/rho_gas
    V_gas = v_gas*m_gas

    #liquid node
    T_liq = T_sat
    m_liq = m_ox-m_gas
    v_liq = 1/rho_liq
    V_liq = v_liq*m_liq

    # NOTE: ASSUME TANK WALL AT T_ATM, MAKE SURE THIS MATCHES REAL WORLD INITIAL CONDITIONS
    T_wall_liq = T_atm
    T_wall_gas = T_atm

    y_0 = np.array([
        T_liq,
        T_gas,
        m_liq,
        m_gas,
        T_wall_liq,
        T_wall_gas,
    ], dtype=np.float64)
    return y_0, rho_liq, rho_gas


def build_run_tank_initial_state(T_atm, P_tank, V_tank, m_liq=1e-2): #m_liq is seed mass

    T_liq = T_sat_table.lookup(P_tank)
    T_gas = T_atm #NOTE: gas solved at P_atm, T_atm, don't know mass need to iter solve


    #NOTE: liq node at quality 0 sat phase
    rho_liq = rho_sat_liq_table.lookup(P_tank)

    v_liq = 1/rho_liq
    V_liq = v_liq*m_liq
            

    m_gas_guess = 0.5 #kg NOTE: initial guess for secant
    V_tank_remaining = V_tank - V_liq

    sol = root_scalar(
        lambda m_gas_guess: P_error(m_gas_guess, V_tank_remaining, T_atm, P_atm),
        method="secant",
        x0=m_gas_guess,
        x1=m_gas_guess * 0.99,
        xtol=1e-9,
        maxiter=100
    )
    if not sol.converged:
        raise RuntimeError("P_dot_error secant solver did not converge")
    
    m_gas = sol.root

    V_gas = V_tank - V_liq
    rho_gas = m_gas/V_gas

    # NOTE: ASSUME TANK WALL AT T_ATM, MAKE SURE THIS MATCHES REAL WORLD INITIAL CONDITIONS
    T_wall_liq = T_atm
    T_wall_gas = T_atm
    

    y_0 = np.array([
        T_liq,
        T_gas,
        m_liq,
        m_gas,
        T_wall_liq,
        T_wall_gas,
    ], dtype=np.float64)
    return y_0, rho_liq, rho_gas

def build_tank_initial_input(P_atm): #same for both!
    
    P_downstream = float(P_atm)
    m_dot_liq_in = 0.0
    h_liq_in = 0.0

    u_0 = np.array([
        P_downstream,
        m_dot_liq_in,
        h_liq_in,
    ], dtype=np.float64)

    return u_0

def update_tank_initial_input(u, P_downstream, m_dot_liq_in, h_liq_in): 
    u[0] = P_downstream
    u[1] = m_dot_liq_in
    u[2] = h_liq_in


def build_tank_param(
        V_tank,
        diam_in,
        diam_out,
        rho_atm,
        T_atm,
        P_atm,
        k_w,
        rho_wall,
        Cd_inj,
        A_inj,
        Cd_vent,
        A_vent,
        CW,
    ):
    
    Tank_Inner_Area, Tank_xsec_Area, height_tank, m_tank = make_tank_geometry_and_consts(diam_in, diam_out, V_tank, rho_wall)

    p = np.array([
        V_tank, 
        diam_in, 
        diam_out, 
        Tank_Inner_Area, 
        height_tank, 
        rho_atm, 
        T_atm, 
        P_atm, 
        k_w, 
        Tank_xsec_Area,
        m_tank, 
        rho_wall, 
        Cd_inj, 
        A_inj, 
        Cd_vent, 
        A_vent, 
        CW, 
    ],dtype=np.float64)
    return p


def build_guess(
        rho_liq_prev, 
        rho_gas_prev, 
        m_dot_inj_out_prev = 0.0, 
        V_dot_liq_prev = 1E-8,
        P_dot_err_tolerance=1.0e-9,
    ):

    guess = np.array([
        rho_liq_prev, 
        rho_gas_prev, 
        m_dot_inj_out_prev,
        V_dot_liq_prev,
        P_dot_err_tolerance,
    ],dtype=np.float64)
    return guess

def update_guess(guess, rho_liq, rho_gas, m_dot_inj_out): # V_dot_liq, P_dot_err_tolerance):
    guess[0] = rho_liq
    guess[1] = rho_gas
    guess[2] = m_dot_inj_out
    #guess[3] = V_dot_liq
    #guess[4] = P_dot_err_tolerance


def resolve_tank_outputs(out, t, y, timestep, u, p, guess):
    

    T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas = y
    P_downstream, m_dot_liq_in, h_liq_in = u

    (
        V_tank,
        diam_in,
        diam_out,
        Tank_Inner_Area,
        height_tank,
        rho_atm,
        T_atm,
        P_atm,
        k_w,
        Tank_xsec_Area,
        m_tank,
        rho_wall,
        Cd_inj,
        A_inj,
        Cd_vent,
        A_vent,
        CW,
    ) = p

    rho_liq_prev, rho_gas_prev, m_dot_inj_out_prev, V_dot_liq_prev, P_dot_err_tolerance = guess

    rho_liq, rho_gas, P_tank, V_liq, V_gas = solve_thermo_params(
        T_liq, T_gas, m_liq, m_gas, rho_liq_prev, rho_gas_prev, V_tank
    )

    h_liq = h_liq_table.lookup(rho_liq, T_liq)

    m_dot_inj_out = -spi_model(Cd_inj, A_inj, P_tank, P_downstream, rho_liq)

    if t <= (timestep * 2.0):
        m_dot_inj_out = 0.5 * m_dot_inj_out
    else:
        m_dot_inj_out = 0.5 * m_dot_inj_out + 0.5 * m_dot_inj_out_prev

    m_dot_vent = spi_model(Cd_vent, A_vent, P_tank, P_downstream, rho_gas)

    out[0]   = P_tank
    out[1]  = rho_liq
    out[2]  = rho_gas
    out[3]    = V_liq
    out[4]    = V_gas
    out[5]    = h_liq
    out[6] = m_dot_inj_out
    out[7] = m_dot_vent

    update_guess(guess, rho_liq, rho_gas, m_dot_inj_out)

    return out





def rhs(t, y, timestep, u, p, guess):
    
    dydt = np.empty(6, dtype=np.float64)

    P_downstream, m_dot_liq_in, h_liq_in = u
    V_tank, diam_in, diam_out, Tank_Inner_Area, height_tank, rho_atm, T_atm, P_atm, k_w, Tank_xsec_Area, m_tank, rho_wall, Cd_inj, A_inj, Cd_vent, A_vent, CW = p
    T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas = y  # Unpack state variables
    rho_liq_prev, rho_gas_prev, m_dot_inj_out_prev, V_dot_liq_prev, P_dot_err_tolerance = guess

    
    rho_liq, rho_gas, P_tank, V_liq, V_gas = solve_thermo_params(T_liq, T_gas, m_liq, m_gas, rho_liq_prev, rho_gas_prev, V_tank)



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
    m_dot_inj = spi_model(Cd_inj, A_inj, P_tank, P_downstream, rho_liq)

    # Also mass transfer from vent
    m_dot_vent = spi_model(Cd_vent, A_vent, P_tank, P_downstream, rho_gas)

    

    # Apply mass flow rate smoothing like equilibrium tank
    if t <= (timestep*2):   # first step
        m_dot_inj = 0.5 * m_dot_inj
        #TODO? m_dot_vent = 0#  0.5 * m_dot_vent
    else:
        m_dot_inj = 0.5 * m_dot_inj + 0.5 * (-m_dot_inj_out_prev)
        #TODO? m_dot_vent = 0.5* m_dot_vent + 0.5 * (-m_dot_vent_prev)




    
    # Heat transfer (2) from saturated surface to gas                       (T_1, T_2, P_tank, rho_2, c, n, tank_diam, fluid)
    # L = tank inner diam , Area of circle x section
    T_film_gas = ((T_sat + T_gas)/2 )
    Q_dot_sat_surf_to_gas = solve_Q_dot_natural_convection_gas(rho_gas, T_sat, T_gas, T_film_gas, P_tank, 0.15, 0.333, diam_in, Tank_Inner_Area, "N2O" ) #relative to gas cv
    
    # Heat transfer (3)  from liq to saturated surface (sat surface assumed to be a liquid with quality 0)
    T_film_liq = ((T_sat + T_liq)/2 )
    Q_dot_liq_to_sat_surf = (E)*solve_Q_dot_natural_convection_liq(rho_liq, T_liq, T_sat, T_film_liq, P_tank, 0.15, 0.333, diam_in, Tank_Inner_Area, "N2O" ) #relative to liq cv
    #NOTE:CORRECTION FACTOR for nitrous oxide heat transfer E = (E) to account for blowing as per [7],[8]

    # Mass transfer (3) by evaporation 
    m_dot_evap = solve_m_dot_evap(h_sat_gas, h_sat_liq, h_liq, Q_dot_liq_to_sat_surf, Q_dot_sat_surf_to_gas)

    # Mass transfer (2) by condensation
    m_dot_cond = solve_m_dot_condensed(T_gas, P_tank, V_gas)


    # Net Mass Transfer of Liquid and Gas CV
    m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj, m_dot_vent, m_dot_liq_in)
    
    #then solve the height of the gas wall
    h_gas_wall = V_gas / Tank_Inner_Area
    V_gas_wall = Tank_xsec_Area*h_gas_wall
    m_gas_wall = rho_wall*V_gas_wall

    h_liq_wall = height_tank - h_gas_wall
    V_liq_wall = Tank_xsec_Area*h_liq_wall #0.25*np.pi*((diam_out**2)-(diam_in**2))*h_liq_wall
    m_liq_wall = rho_wall*V_liq_wall


    # Heat transfer (4) [natural convection] from liq wall to liq
    Q_dot_liq_wall_to_liq = solve_Q_dot_natural_convection_liq(rho_liq, T_wall_liq, T_liq, T_liq, P_tank, 0.021, 0.4, h_liq_wall, (np.pi*diam_in*h_liq_wall + Tank_Inner_Area), "N2O" ) #relative to liq cv       
    # Heat transfer (5) [natural convection] from gas to gas wall
    Q_dot_gas_wall_to_gas = solve_Q_dot_natural_convection_gas(rho_gas, T_wall_gas, T_gas, T_gas, P_tank, 0.021, 0.4, h_gas_wall, (np.pi*diam_in*h_gas_wall + Tank_Inner_Area), "N2O" ) #relative to gas cv
    

    Q_dot_liq = Q_dot_liq_wall_to_liq - Q_dot_liq_to_sat_surf
    Q_dot_gas = Q_dot_gas_wall_to_gas + Q_dot_sat_surf_to_gas


    print("Q_dot_gas: ", Q_dot_gas, Q_dot_gas_wall_to_gas, Q_dot_sat_surf_to_gas)


    #NOTE: USE ambient properties for air, T_2 will be respective wall temperature (RK var)
    # (6) [natural convection] from atm to liq wall
    Q_dot_atm_to_liq_wall = solve_Q_dot_natural_convection_gas(rho_atm, T_atm, T_wall_liq, T_atm, P_atm, 0.59, 0.25, height_tank, (np.pi*diam_out*h_liq_wall + Tank_Inner_Area), "Air") #relative to wall_liq cv
    # (7) [natural convection] from atm to gas wall
    Q_dot_atm_to_gas_wall = solve_Q_dot_natural_convection_gas(rho_atm, T_atm, T_wall_gas, T_atm, P_atm, 0.59, 0.25, height_tank, (np.pi*diam_out*h_gas_wall + Tank_Inner_Area), "Air") #relative to wall_gas cv
    
    # (8) [conduction] from liq wall to gas wall 
    Q_dot_liq_wall_to_gas_wall = solve_Q_dot_conduction( (T_wall_liq-T_wall_gas), height_tank, k_w, Tank_xsec_Area) #relative to wall_liq cv

    # Iteratively solve change in CV Volume
    V_dot_liq = V_dot_liq_prev #initial guess for dV_dt_liq
    


    
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
        xtol=P_dot_err_tolerance,
        maxiter=100
    )

    if not sol.converged:
        raise RuntimeError("P_dot_error secant solver did not converge")

    V_dot_liq = sol.root        



    ### Wall nodes:
    
    height_dot = V_dot_liq / Tank_Inner_Area

    m_dot_liq_wall = rho_wall * Tank_xsec_Area * height_dot
    m_dot_gas_wall = -m_dot_liq_wall


    #NOTE: for T_dot_wall_liq:  IN: (6)      OUT: (4) AND (8)
    T_dot_wall_liq = ( Q_dot_atm_to_liq_wall - Q_dot_liq_wall_to_liq - Q_dot_liq_wall_to_gas_wall + m_dot_liq_wall*CW*(T_wall_liq - T_wall_gas) ) / (CW*m_liq_wall)

    #NOTE: for T_dot_wall_vap:  IN: (7) and (8)      OUT: (5)
    T_dot_wall_gas = ( Q_dot_atm_to_gas_wall - Q_dot_gas_wall_to_gas + Q_dot_liq_wall_to_gas_wall + m_dot_gas_wall*CW*(T_wall_liq - T_wall_gas) ) / (CW*m_gas_wall)

    T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_liq, P_tank, m_liq, m_gas, V_liq, V_gas, cv_liq, cv_gas, h_gas, h_sat_gas, h_sat_liq, h_liq, u_liq, u_gas, du_drho_const_T_liq, du_drho_const_T_gas, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, m_dot_vent, m_dot_liq_in, h_liq_in, True)




    dydt[0] = T_dot_liq
    dydt[1] = T_dot_gas
    dydt[2] = m_dot_liq
    dydt[3] = m_dot_gas
    dydt[4] = T_dot_wall_liq
    dydt[5] = T_dot_wall_gas
    return dydt







def main():
    t = 0

    ### NITROUS BOTTLE ###

    botl_P_tank = 45e5 #Pa #TODO: check this across envelope, worst case fill time at T_min

    botl_m_nos = 29.4835 #kg #NOTE: 65LB converted to kg
    botl_V_tank = 0.1 #m^3 #NOTE: PLACEHOLDER VAL, SHOULD BE: this is for a Q type, can either buy K type or Q type from Linde

    botl_diam_in = (0.0254*7) #m #from Linde website
    botl_diam_out = botl_diam_in + 2*(0.0254*0.25) #m #guessing, as stated before wall node

    # i believe bottle is al but either way i dont believe wall nodes at our scale matter much
    botl_rho_wall = 2770 #kg/m^3
    botl_k_w = 237 #W/(m K)

    botl_Cd_inj = 0.4
    botl_A_inj = 0.25*np.pi*(0.0254*(.2))**2

    #no venting, placeholder
    botl_Cd_vent = 0
    botl_A_vent = 0

    
    botl_y, botl_rho_liq, botl_rho_gas = build_bottle_initial_state(botl_P_tank, botl_V_tank, botl_m_nos, T_atm)
    botl_p = build_tank_param(botl_V_tank, botl_diam_in, botl_diam_out, rho_atm, T_atm, P_atm, botl_k_w,botl_rho_wall, botl_Cd_inj, botl_A_inj, botl_Cd_vent, botl_A_vent, CW)
    botl_u = build_tank_initial_input(P_atm)
    botl_guess = build_guess(botl_rho_liq, botl_rho_gas)
    botl_out = np.empty(8, dtype=np.float64)




    ### S1 RUN TANK ### - GND % fill
    full_pcnt_fill = 0.90
    #TODO: UPDATE

    run_tank_m_nos = 0 #kg NOTE: placeholder, solved to obtain density for thermo state at P_atm, T_atm 
    run_tank_P_tank = P_atm #Pa #NOTE: starts as empty, assuming initial condition of tank is filled w nitrous @ P_atm
    run_tank_V_tank = 0.0208 #m^3

    run_tank_diam_in = 0.1282 #m 
    run_tank_diam_out = 0.1413 #m 
    run_tank_rho_wall = 2770 #kg/m^3
    run_tank_k_w = 237 #W/(m K)

    run_tank_Cd_vent = 0.65 #using this from RPE 9th ed p279 
    run_tank_A_vent = 0.25*np.pi*(0.0254*(.002))**2 #NOTE: assuming 

    #no draining, placeholder
    run_tank_Cd_inj = 0.65 #using this from RPE 9th ed p279 for sharp edge orifice
    run_tank_A_inj = 0.25*np.pi*(0.0254*(.010))**2 #NOTE: assuming 2'" diam outlet

    """
    ### S1 RUN TANK ### - FLIGHT % fill
    full_pcnt_fill = 0.9

    run_tank_m_nos = 13.25 #kg
    run_tank_P_tank = P_atm #Pa #NOTE: starts as empty, assuming initial condition of tank is filled w nitrous @ P_atm
    run_tank_V_tank = 0.0208 #m^3

    run_tank_diam_in = 0.1282 #m 
    run_tank_diam_out = 0.1413 #m 
    run_tank_rho_wall = 2770 #kg/m^3
    run_tank_k_w = 237 #W/(m K)

    run_tank_Cd_vent = 0.65 #using this from RPE 9th ed p279 
    run_tank_A_vent = 0.25*np.pi*(0.0254*(.002))**2 #NOTE: assuming 

    #no draining, placeholder
    run_tank_Cd_inj = 0.57
    run_tank_A_inj = 0.25 * 3.14159 * (42) * (1.5e-3)**2 
    """


    run_tank_y, run_tank_rho_liq, run_tank_rho_gas = build_bottle_initial_state(run_tank_P_tank, run_tank_V_tank, 1, T_atm)#build_run_tank_initial_state(T_atm, run_tank_P_tank, run_tank_V_tank)
    run_tank_p = build_tank_param(run_tank_V_tank, run_tank_diam_in, run_tank_diam_out, rho_atm, T_atm, P_atm, run_tank_k_w, run_tank_rho_wall, run_tank_Cd_inj, run_tank_A_inj, run_tank_Cd_vent, run_tank_A_vent, CW)
    run_tank_u = build_tank_initial_input(P_atm)
    run_tank_guess = build_guess(run_tank_rho_liq, run_tank_rho_gas)
    run_tank_out = np.empty(8, dtype=np.float64)

    
    
    try:
        start_time = time.time()  # Start timer

        m_vented = 0.0

        pcnt_fill = 0.0

        log_every = 60
        step_idx = 0

        while( pcnt_fill < full_pcnt_fill): #t <=1000*TIMESTEP): #pcnt_fill < full_pcnt_fill): #t <= 3.0):
            
            print("botle")
            botl_y = rk4_step(rhs, t, botl_y, TIMESTEP, botl_u, botl_p, botl_guess)

            resolve_tank_outputs(botl_out, t, botl_y, TIMESTEP, botl_u, botl_p, botl_guess)
            #m_dot_inj_out, h_liq = botl.inst(P_run_tank, 0, 0)

            update_tank_initial_input(run_tank_u, 0.0, botl_out[6], botl_out[5])

            #add rocketprops check here

            print("run tank")
            run_tank_y = rk4_step(rhs, t, run_tank_y, TIMESTEP, run_tank_u, run_tank_p, run_tank_guess)

            resolve_tank_outputs(run_tank_out, t, run_tank_y, TIMESTEP, run_tank_u, run_tank_p, run_tank_guess)

            update_tank_initial_input(botl_u, run_tank_out[0], 0.0, 0.0)


            ###sol vented mass in that step
            m_vented += TIMESTEP * run_tank_out[7]


            #march timestep
            t+=TIMESTEP 
            

            #recalc pcnt_fill and P_run_tank for next timestep
            pcnt_fill = run_tank_out[3]/run_tank_p[0]

            # logging

            if step_idx % log_every == 0:
                append_dat_resolved(t, botl_dat, botl_y, botl_out)
                append_dat_resolved(t, run_tank_dat, run_tank_y, run_tank_out)

                #print(f"m_dot_into_run_tank: {botl_out[6]}")
                print(f"%fill: {pcnt_fill:.2f}, m_vented: {m_vented:.5f}")

                #print(f"in bottle at t = {t}, final y: ", botl.T_liq, botl.T_gas, botl.m_liq, botl.m_gas, botl.T_wall_liq, botl.T_wall_gas, "\n")
                #print(f"in run tank at t = {t}, final y: ", run_tank.T_liq, run_tank.T_gas, run_tank.m_liq, run_tank.m_gas, run_tank.T_wall_liq, run_tank.T_wall_gas, "\n")

            step_idx += 1


            

        end_time = time.time()  # End timer
        #print(f"\nTotal real time for simulation: {end_time - start_time:.3f} seconds")


        m_bottle_remaining = botl_y[2] + botl_y[3]

        #my bad figma moment
        print("\n\n\n\n\n****************************************************")
        print("*          SOLUTION TERMINATED NORMALLY            *")
        print("****************************************************")
        print(f"*  REAL SIMULATION TIME [s]     = {end_time - start_time:12.5f}    *")
        print(f"*  FILL TIME [s]                = {t:12.5f}    *")
        print(f"*  VENTED MASS [kg]             = {m_vented:12.5f}    *")
        print(f"*  BOTTLE MASS REMAINING [kg]   = {m_bottle_remaining:12.5f}    *")
        print("****************************************************")

    except Exception as e:
        traceback.print_exc()

    plot_dat(botl_dat, "BOTTLE")
    plot_dat(run_tank_dat, "RUN TANK")
    plt.show()

    pass



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

def append_dat_resolved(t, dat, y, out):
    T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas = y

    dat.time_arr.append(t)
    dat.P_tank_arr.append(out[0])
    dat.m_liq_arr.append(m_liq)
    dat.m_gas_arr.append(m_gas)

    dat.T_liq_arr.append(T_liq)
    dat.T_gas_arr.append(T_gas)

    try:
        dat.T_sat_arr.append(T_sat_table.lookup(out[0]))
    except Exception:
        dat.T_sat_arr.append(0.0)

    dat.T_liq_wall_arr.append(T_wall_liq)
    dat.T_gas_wall_arr.append(T_wall_gas)


def plot_dat(dat, fig_name):
    plt.figure()

    plt.suptitle(fig_name, fontsize=16)

    plt.subplot(1,3,1)
    plt.plot(dat.time_arr,dat.P_tank_arr,label = "model_tank", color = "b" )

    #plt.plot(published_time_arr,published_P_tank_arr,label = "published_tank", color = "dodgerblue" )
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (Pa)')
    plt.title('Pressure vs. Time')
    plt.legend()
    plt.grid(True)

    plt.subplot(1,3,2)
    plt.plot(dat.time_arr,dat.m_liq_arr, label = "model_liquid", color = "b" )
    plt.plot(dat.time_arr,dat.m_gas_arr, label = "model_gas", color = "darkorange" )


    plt.xlabel('Time (s)')
    plt.ylabel('Mass (kg)')
    plt.title('Mass vs. Time')
    plt.legend()
    plt.grid(True)

    plt.subplot(1,3,3)
    plt.plot(dat.time_arr,dat.T_liq_arr, label = "model_liquid", color = "b" )
    plt.plot(dat.time_arr,dat.T_gas_arr, label = "model_gas", color = "darkorange" )
    plt.plot(dat.time_arr,dat.T_sat_arr, label = "model_T_sat", color =  "g" )


    plt.plot(dat.time_arr,dat.T_liq_wall_arr, label = "model_WALL liquid", color = "r" )
    plt.plot(dat.time_arr,dat.T_gas_wall_arr, label = "model_WALL gas", color = "mediumorchid" )
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature vs. Time')
    plt.legend()
    plt.grid(True)


if __name__ == "__main__":
    #profiler = cProfile.Profile()
    #profiler.enable()
    main()
    #profiler.disable()

    #stats = pstats.Stats(profiler)
    #stats.strip_dirs()
    #stats.sort_stats("cumtime")   # or "tottime"
    #for func, stat in stats.stats.items():
    #    cc, nc, tt, ct, callers = stat
    #    print(f"{func}  ncalls={nc}  tottime={tt:.9f}  cumtime={ct:.9f}")
    #stats.print_stats(60)         # top n functions