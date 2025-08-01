from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo import Chemical
from thermo.eos import PR
import matplotlib.pyplot as plt
import traceback
import time

#from ..thermo_property_lookup.get_n2o_viscosity import get_n2o_viscosity
#python3 -m src.models.Karabeyoglu-Zilliac-Zimmerman-nos-tank

### this is to test:
LOOKUP_TIME = 0

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

H_OFFSET = 0#4.56e4

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



### ig polynomials here!!!!
def gas_dynamic_visc_polynomial(T):
    # Polynomial coefficients
    A = 2.1150E-6
    B = 0.46420
    C = 305.70
    D = 0.0

    # Apply temperature limits
    if 182 < T and T < 1000:
        dvisc = A*T**B / ( 1 + C/T + D/T**2)
        return dvisc
    raise ValueError("Temperature outside of function bounds!")

def liq_dynamic_visc_polynomial(T):
    # Polynomial coefficients
    A = 0.001877085
    B = -1.1864E-5  #NOTE: Thesis reported B = 1.1864E-5 but this will return a dynamic viscosity orders of magnitude higher than expected. From testing this function found it was likely signed wrong in the thesis.
    C = 1.928E-8 

    # Apply temperature limits
    if 182 < T and T < 1000:
        dvisc = A + B*T + C*T**2
        return dvisc
    raise ValueError("Temperature outside of function bounds!")


def solve_cp_ig_polynomial(T):
    # Polynomial coefficients
    A = 21.62
    B = 72.81
    C = -57.78  
    D = 18.3
    E = 0.0

    # Apply temperature limits
    if 150 < T and T < 310:
        T_reduced = T / 1000
        cp_ig = (A + B * T_reduced + C * T_reduced**2 + D * T_reduced**3 + E / (T_reduced**2)) / MW  # J/(kg K)
        return cp_ig
    raise ValueError("Temperature outside of function bounds!")

def solve_cv_ig_polynomial(T):
    cv_ig = solve_cp_ig_polynomial(T) - R_U
    return cv_ig

def analytical_integration_ig_enthalpy(T_REF, T):
    # Polynomial coefficients
    A = 21.62
    B = 72.81
    C = -57.78  
    D = 18.3
    E = 0.0

    # Apply temperature limits
    if 150 < T and T < 310:
        h_ig = (T*(12000000000000000*E - T_REF**2*(12000000000*A + 6000000*B*T_REF + 4000*C*T_REF**2 + 3*D*T_REF**3)) + T_REF*(-12000000000000000*E + T**2*(12000000000*A + 6000000*B*T + 4000*C*T**2 + 3*D*T**3)))/(12000000000*MW*T*T_REF)
        return h_ig
    raise ValueError("Temperature outside of function bounds!")


def analytical_integration_ig_int_energy(T_REF, T):
    u_ig = analytical_integration_ig_enthalpy(T_REF, T) - R_U*T
    return u_ig


def solve_du_drho_const_T_liq_central_diff(T_liq, rho_liq, delta_rho):
    # my bad guys its numerical bc kept messing up the calculus
    rho_0 = rho_liq - delta_rho
    vm_0 = MW/rho_0 
    preos_l_0 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, V= vm_0)
    u_dep_0 = preos_l_0.U_dep_l

    rho_1 = rho_liq + delta_rho
    vm_1 = MW/rho_1
    preos_l_1 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, V= vm_1)
    u_dep_1 = preos_l_1.U_dep_l

    central_difference_du_drho_const_T_liq =  ( (u_dep_1/MW) - (u_dep_0/MW) ) / (rho_1-rho_0) #since u_ig = f(T) only, the ideal gas component is the same so take difference in departure function
    return central_difference_du_drho_const_T_liq


def solve_du_drho_const_T_gas_central_diff(T_gas, rho_gas, delta_rho):
    # my bad guys its numerical bc kept messing up the calculus
    rho_0 = rho_gas - delta_rho
    vm_0 = MW/rho_0 
    preos_g_0 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, V= vm_0)
    u_dep_0 = preos_g_0.U_dep_g

    rho_1 = rho_gas + delta_rho
    vm_1 = MW/rho_1
    preos_g_1 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, V= vm_1)
    u_dep_1 = preos_g_1.U_dep_g

    central_difference_du_drho_const_T_gas =  ( (u_dep_1/MW) - (u_dep_0/MW) ) / (rho_1-rho_0) #since u_ig = f(T) only, the ideal gas component is the same so take difference in departure function
    return central_difference_du_drho_const_T_gas



#NOTE: SIGN CONVENTION: Q_dot dir: (+) T_1 --> T_2 (f is fluid)
def solve_Q_dot_natural_convection_liq(T_1, T_2, T_f, P_f, rho_f, c, n, L, Area, fluid):

    if fluid == "N2O":
        n2o = Chemical('N2O', T=T_f, P=P_f) 
        k_f = n2o.kl # Conductivity W/(m K)
        dyn_visc_f = liq_dynamic_visc_polynomial(T_f)
        visc_f = dyn_visc_f/rho_f
        #visc_f = get_n2o_viscosity(T_f, P_f, "liquid") # Kinematic viscosity (m^2/s)

        cp_ig = solve_cp_ig_polynomial(T_f)

        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_f, P=P_f)
        cp_f = (preos_l.Cp_dep_l/MW + cp_ig) #J/K 

        dV_dT_P = preos_l.dV_dT_l
        beta = (1/preos_l.V_l)*dV_dT_P

    elif fluid == "Air":
        k_f = CP.PropsSI('L', 'T', T_f, 'P', P_f, 'Air')  # Conductivity W/(m K)
        dyn_visc_f = CP.PropsSI('V', 'T', T_f, 'P', P_f, 'Air')  # Dynamic viscosity (Pa s)
        visc_f = dyn_visc_f/rho_f
        
        cp_f = CP.PropsSI("Cpmass", "T", T_f, "P", P_f, "Air")
        beta = CP.PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T_f, "P", P_f, "Air")

    Gr = ((L**3)*g*beta*np.abs(T_2 - T_1) ) / (visc_f**2)
    Pr = (cp_f*visc_f)/ k_f
    X = Gr*Pr

    h = c * (k_f/L) * X**n

    Q_dot = h*Area*(T_1-T_2)

    return Q_dot #NOTE: Q_dot + going into (2)



#NOTE: SIGN CONVENTION: Q_dot dir: (+) T_1 --> T_2 
def solve_Q_dot_natural_convection_gas(T_1, T_2, T_f, P_f, rho_f, c, n, L, Area, fluid): #BUG: potential mistake, solving _1 properties with _2 inputs, double check this is likely a mistake

    if fluid == "N2O":
        n2o = Chemical('N2O', T=T_f, P=P_f)  #TODO: units here!!!
        k_f = n2o.kg # Conductivity W/(m K)
        dyn_visc_f = gas_dynamic_visc_polynomial(T_f)
        visc_f = dyn_visc_f/rho_f
        #visc_f = get_n2o_viscosity(T_f, P_f, "vapor") # Kinematic viscosity (m^2/s)

        cp_ig = solve_cp_ig_polynomial(T_f)

        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_f, P=P_f)
        cp_f = (preos_g.Cp_dep_g/MW + cp_ig)
        
        dV_dT_P = preos_g.dV_dT_g 
        beta = (1/preos_g.V_g)*dV_dT_P

    elif fluid == "Air":
        k_f = CP.PropsSI('L', 'T', T_f, 'P', P_f, 'Air')  # Conductivity W/(m K)
        dyn_visc_f = CP.PropsSI('V', 'T', T_f, 'P', P_f, 'Air')  # Dynamic viscosity (Pa s)
        visc_f = dyn_visc_f/rho_f 
        
        cp_f = CP.PropsSI("Cpmass", "T", T_f, "P", P_atm, "Air")
        beta = CP.PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T_f, "P", P_f, "Air")

    Gr = ((L**3)*g*beta*np.abs(T_2 - T_1) ) / (visc_f**2)
    Pr = (cp_f*visc_f)/ k_f
    X = Gr*Pr

    h = c * (k_f/L) * X**n

    Q_dot = h*(Area)*(T_1-T_2)

    return Q_dot



def solve_Q_dot_conduction(delta_T, h_tank, k_w, diam_in, diam_out):

    L_w_cond = 0.5*h_tank
    Q_dot_conduction = k_w *(delta_T)*(0.25*np.pi*((diam_out**2)-(diam_in**2)))/L_w_cond
    
    return Q_dot_conduction



#TODO: update with other injector model once we get this thing up
def spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit):

    """
    #for tomacz test case!!!
    pipe_inj_time = [ 0, 0.25, 1.1, 1.5,4]
    pipe_inj_m_dot = [ (-50/1000), (-43/1000), (-41.8/1000), (-36/1000), (-22/1000)]

    m_dot_spi = np.interp(LOOKUP_TIME, pipe_inj_time , pipe_inj_m_dot)
    """
#NOTE: PIPING M_DOT TO ISOLATE REST OF MODEL FOR DEBUG
    m_dot_spi = -3.75
    return m_dot_spi


def solve_m_dot_evap( T_gas, T_liq, P_tank, Q_dot_liq_to_sat_surf, Q_dot_sat_surf_to_gas):
    m_dot_evap = 0
    #print("evap heat transfer rates: ", (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas), Q_dot_liq_to_sat_surf, Q_dot_sat_surf_to_gas)
    if (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) > 0:

        #old, this is wrong and does not match the thesis!
        """
        h_ig_liq = analytical_integration_ig_enthalpy(T_REF, T_liq)

        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
        h_liq = preos_l.H_dep_l/MW + h_ig_liq 

        h_ig_gas = analytical_integration_ig_enthalpy(T_REF, T_gas)

        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
        h_gas = preos_g.H_dep_g/MW + h_ig_gas 

        m_dot_evap = (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) / (h_gas-h_liq) 
        """

        h_ig_liq = analytical_integration_ig_enthalpy(T_REF, T_liq)

        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
        h_liq = preos_l.H_dep_l/MW + h_ig_liq + H_OFFSET

        T_sat = preos_l.Tsat(P_tank)

        h_ig_sat = analytical_integration_ig_enthalpy(T_REF, T_sat)

        preos_sat = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_sat, P=P_tank)
        h_sat_liq = preos_sat.H_dep_l/MW + h_ig_sat + H_OFFSET 
        h_sat_gas = preos_sat.H_dep_g/MW + h_ig_sat + H_OFFSET

        m_dot_evap = (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) / ( (h_sat_gas-h_sat_liq) + (h_sat_liq - h_liq)  )  #if this doesnt work check functions and if they are as reusable as i treat them

        #print("sign convention in m_dot_evap Q: ", Q_dot_liq_to_sat_surf, - Q_dot_sat_surf_to_gas)

    return m_dot_evap


def solve_m_dot_condensed(T_gas, T_liq, P_tank, V_gas, t):
    m_dot_cond = 0
    
    preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
    P_sat_g = preos_g.Psat(T_gas)

    if (P_tank > P_sat_g):
        m_dot_cond = ((P_tank-P_sat_g)*V_gas*MW)/( preos_g.Z_g*(R_U/MW)*T_gas*(TIMESTEP) )  

#NOTE: CONDENSATION FROM PREOS
    #if p_tank > p_sat_gas, then condensation to enforce equilibrium

###NOTE: try m_dot_cond = 0
    m_dot_cond = 0
    return m_dot_cond



def solve_U_dot_liq(T_liq, T_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq , Q_dot_net, h_liq_prev, h_gas_prev):

    
    preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    T_sat = preos_l.Tsat(P_tank)
    preos_sat = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_sat, P=P_tank)
    #preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)

    h_ig_liq = analytical_integration_ig_enthalpy(T_REF, T_liq)
    h_liq = preos_l.H_dep_l/MW + h_ig_liq + H_OFFSET

    h_ig_sat = analytical_integration_ig_enthalpy(T_REF, T_sat)
    h_sat_gas = preos_sat.H_dep_g/MW + h_ig_sat + H_OFFSET
    h_sat_liq = preos_sat.H_dep_l/MW + h_ig_sat + H_OFFSET

    """
    h_ig_gas = analytical_integration_ig_enthalpy(T_REF, T_gas)
    h_gas = preos_g.H_dep_g/MW + h_ig_gas + H_OFFSET
    """

    # test U_dot_liq = -m_dot_inj*( (h_liq-h_liq_prev ) ) - m_dot_evap*( (h_gas-h_gas_prev) ) -P_tank*V_dot_liq + Q_dot_net
    #U_dot_liq = m_dot_inj*h_liq  + m_dot_cond*( 0 ) -P_tank*V_dot_liq + Q_dot_net

    U_dot_liq = m_dot_inj*(h_liq*.08) - m_dot_evap*(h_sat_gas - h_sat_liq) + m_dot_cond*( 0 ) -P_tank*V_dot_liq + Q_dot_net
    #NOTE: m_dot_inj term might be a spatial difference in enthalpy, keep this in mind if more issues

    #print("U_dot_liq term signs: ",U_dot_liq, m_dot_inj*h_liq, - m_dot_evap*(h_sat_gas - h_sat_liq), m_dot_cond*( 0 ), -P_tank*V_dot_liq , Q_dot_net)

    return U_dot_liq



def solve_U_dot_gas(T_liq, T_gas, P_tank, m_dot_evap, m_dot_cond, V_dot_gas, Q_dot_net, h_liq_prev, h_gas_prev):

    preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    T_sat = preos_l.Tsat(P_tank)
    preos_sat = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_sat, P=P_tank)
    preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)

    """
    h_ig_liq = analytical_integration_ig_enthalpy(T_REF, T_liq)
    h_liq = preos_l.H_dep_l/MW + h_ig_liq + H_OFFSET
    """

    h_ig_sat = analytical_integration_ig_enthalpy(T_REF, T_sat)
    h_sat_gas = preos_sat.H_dep_g/MW + h_ig_sat + H_OFFSET
    #h_sat_liq = preos_sat.H_dep_l/MW + h_ig_sat + H_OFFSET

    h_ig_gas = analytical_integration_ig_enthalpy(T_REF, T_gas)
    #h_gas = preos_g.H_dep_g/MW + h_ig_gas + H_OFFSET

    U_dot_gas = m_dot_evap*h_sat_gas - m_dot_cond*( (0) ) - P_tank*V_dot_gas + Q_dot_net 
#BUG: h_sat_gas has a negative sign which is messing up U_dot_gas
    #print("U_dot_gas term signs: ", U_dot_gas, m_dot_evap*h_sat_gas, - m_dot_cond*( (0) ), - P_tank*V_dot_gas , Q_dot_net )
    #print( m_dot_evap, U_dot_gas, h_sat_gas)

    return U_dot_gas



def solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj):
    m_dot_liq = -m_dot_evap + m_dot_cond + m_dot_inj
    m_dot_gas = m_dot_evap - m_dot_cond #convert sign convention from liq cv to gas cv

    return m_dot_liq, m_dot_gas



def V_tank_error(P_guess, T_liq, T_gas, m_liq, m_gas, V_tank):

    preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_guess)
    rho_liq = preos_l.rho_l*MW
    
    preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_guess)
    rho_gas = preos_g.rho_g*MW

    V_tank_est = (m_liq/rho_liq) + (m_gas/rho_gas)

    return V_tank_est - V_tank 



def solve_thermo_params(T_liq, T_gas, m_liq, m_gas, P_tank_prev, V_tank, volume_err_tolerance):

    P_tank = P_tank_prev #initial guess for pressure

    while np.abs(V_tank_error(P_tank, T_liq, T_gas, m_liq, m_gas, V_tank) ) > volume_err_tolerance:
        P_tank = secant((lambda P: V_tank_error(P, T_liq, T_gas, m_liq, m_gas, V_tank)), P_tank)

    preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    rho_liq = preos_l.rho_l*MW
    
    preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
    rho_gas = preos_g.rho_g*MW

    return rho_liq, rho_gas, P_tank



def single_solve_T_dot_liq_gas(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, h_liq_prev, h_gas_prev, u_liq_prev, u_gas_prev, debug_mode):

    m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj)

    V_dot_gas = -V_dot_liq

    d_rho_dt_liq = (1/V_liq)*m_dot_liq -(m_liq/(V_liq**2))*V_dot_liq
    d_rho_dt_gas = (1/V_gas)*m_dot_gas -(m_gas/(V_gas**2))*V_dot_gas

    U_dot_liq = solve_U_dot_liq(T_liq, T_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq, Q_dot_liq, h_liq_prev, h_gas_prev)
    U_dot_gas = solve_U_dot_gas(T_liq, T_gas, P_tank, m_dot_evap, m_dot_cond, V_dot_gas, Q_dot_gas, h_liq_prev, h_gas_prev)
    
    preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
    partial_du_d_rho_const_T_gas = solve_du_drho_const_T_gas_central_diff(T_gas,rho_gas,0.00001)
    
    cv_ig_gas = solve_cv_ig_polynomial(T_gas)
    cv_gas = preos_g.Cv_dep_g/MW + cv_ig_gas

    u_ig_gas = analytical_integration_ig_int_energy(T_REF, T_gas)
    u_gas = preos_g.U_dep_g/MW + u_ig_gas 


    preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    partial_du_d_rho_const_T_liq = solve_du_drho_const_T_liq_central_diff(T_liq,rho_liq,0.00001)

    cv_ig_liq = solve_cv_ig_polynomial(T_liq)
    cv_liq = (preos_l.Cv_dep_l/MW + cv_ig_liq)

    u_ig_liq = analytical_integration_ig_int_energy(T_REF, T_liq)
    u_liq = preos_l.U_dep_l/MW + u_ig_liq 


    T_dot_liq = (1/cv_liq)*( (1/m_liq) * (U_dot_liq - (u_liq*m_dot_liq)) - (partial_du_d_rho_const_T_liq* d_rho_dt_liq) )
    T_dot_gas = (1/cv_gas)*( (1/m_gas) * (U_dot_gas - (u_gas*m_dot_gas)) - (partial_du_d_rho_const_T_gas* d_rho_dt_gas) )


    if debug_mode == True:

        ### check energy balance:
        U_dot_liq_check = m_liq*cv_liq*T_dot_liq + m_liq*(partial_du_d_rho_const_T_liq* d_rho_dt_liq) + u_liq*m_dot_liq
        U_dot_gas_check = m_gas*cv_gas*T_dot_gas + m_gas*(partial_du_d_rho_const_T_gas* d_rho_dt_gas) + u_gas*m_dot_gas


        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
        h_ig_liq = analytical_integration_ig_enthalpy(T_REF, T_liq)
        h_liq = preos_l.H_dep_l/MW + h_ig_liq + H_OFFSET


        print("cons energy CHECK IN T_DOT_EQN: ", (U_dot_liq_check + U_dot_gas_check), " = ", m_dot_inj*h_liq, "expecting 0: ",  (U_dot_liq_check + U_dot_gas_check)-(m_dot_inj*h_liq), (P_tank*(V_dot_liq+(-V_dot_liq))) )


        a = 1
        #print("T_dot_liq: ", T_dot_liq, (U_dot_liq - (u_liq-u_liq_prev)*np.abs(m_dot_evap)) , - (partial_du_d_rho_const_T_liq* d_rho_dt_liq)  )
        #print("U_dot liq/gas: ", U_dot_liq, U_dot_gas)
        #print("looking at T_dot: ", m_dot_evap, T_dot_liq, T_dot_gas)
#TODO: delete at some point

    return T_dot_liq, T_dot_gas 

def P_dot_error(V_dot_guess, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, h_liq_prev, h_gas_prev, u_liq_prev, u_gas_prev):    

    V_dot_gas = -V_dot_guess #guessing for liquid

    m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj) #or put this outside and pass in?


    d_rho_dt_liq = (1/V_liq)*m_dot_liq - (m_liq/(V_liq**2))*V_dot_guess
    d_rho_dt_gas = (1/V_gas)*m_dot_gas - (m_gas/(V_gas**2))*V_dot_gas

    T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_guess, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, h_liq_prev, h_gas_prev, u_liq_prev, u_gas_prev, False)


    preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    partial_dP_dT_const_rho_liq = preos_l.dP_dT_l #- (MW/(rho_liq**2))*preos_l.dP_dV_l #NOTE: deriv at const Vm != deriv const rho ???

    partial_dP_drho_const_T_liq = (preos_l.dP_dV_l)*(-MW/(rho_liq**2)) #converting deriv wrt Vm --> wrt rho

    P_dot_liq = partial_dP_dT_const_rho_liq*T_dot_liq + partial_dP_drho_const_T_liq*d_rho_dt_liq


    preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
    partial_dP_dT_const_rho_gas = preos_g.dP_dT_g #- (MW/(rho_gas**2))*preos_g.dP_dV_g #NOTE: deriv const Vm != deriv const rho ???

    partial_dP_drho_const_T_gas = (preos_g.dP_dV_g)*(-MW/(rho_gas**2))  #converting deriv wrt Vm --> wrt rho

    P_dot_gas = partial_dP_dT_const_rho_gas*T_dot_gas + partial_dP_drho_const_T_gas*d_rho_dt_gas


    #print("P_dot_inside (take last) :", P_dot_liq, P_dot_gas )
    return P_dot_liq - P_dot_gas


class model():

    def __init__(self, timestep, m_nos, P_tank, P_cc, P_atm, T_atm, rho_atm, V_tank, diam_out, diam_in, rho_wall, k_w, Cd_1, A_inj_1, volume_err_tolerance, P_dot_err_tolerance):
        
        self.TIMESTEP = timestep
        self.volume_err_tolerance = volume_err_tolerance
        self.P_dot_err_tolerance = P_dot_err_tolerance

        self.timestep = timestep

        self.T_atm = T_atm
        self.P_atm = P_atm
        self.rho_atm = rho_atm

        #NOTE: Tank geometry and setup
        self.rho_wall = rho_wall
        self.k_w = k_w
        self.diam_out = diam_out
        self.diam_in = diam_in

        #NOTE: System starts with valve closed assuming thermal equillibrium:
        self.m_dot_liq = 0
        self.m_dot_gas = 0

        #injector input
        self.Cd_1 = Cd_1
        self.P_cc = P_cc
        self.A_inj_1 = A_inj_1

        self.P_tank = P_tank

        #### Setup Thermo Properties - Assuming Tank Starts at thermal equilibrium --> Sat Conditions

        self.V_tank = V_tank
        self.m_nos = m_nos

        preos = PR(Tc=TC, Pc=PC, omega=OMEGA, T=self.T_atm, P=self.P_tank)
        rho_sat_gas = (1/preos.V_g_sat(self.T_atm))*MW
        rho_sat_liq = (1/preos.V_l_sat(self.T_atm))*MW
        
        rho_sat_tank = self.m_nos/V_tank

        x_tank = ( (1/rho_sat_tank)-(1/rho_sat_liq) ) / ( (1/rho_sat_gas)-(1/rho_sat_liq) )

        #gas cv
        self.m_gas = x_tank*self.m_nos

        self.T_gas = preos.Tsat(self.P_tank) #- 0.0005 #Perturb to start close to equillibrium but as a gas
        self.T_wall_gas = self.T_gas

        preos = PR(Tc=TC, Pc=PC, omega=OMEGA, T=self.T_gas, P=self.P_tank)
        self.rho_gas = preos.rho_g*MW
        v_gas = 1/self.rho_gas
        self.V_gas = v_gas*self.m_gas

        #liquid cv
        self.m_liq = self.m_nos-self.m_gas
        self.T_liq = preos.Tsat(self.P_tank) #+ 0.0005 #Perturb to start close to equillibrium
        self.T_wall_liq = self.T_liq

        #solve rho_liq!
        self.rho_liq = preos.rho_l*MW
        v_liq = 1/self.rho_liq
        self.V_liq = v_liq*self.m_liq


        self.V_tank = self.V_liq+self.V_gas # "what are you going to do if the aluminum is too small? water it? give it sunlight? let it grow?"
        self.height_tank = self.V_tank/(0.25*np.pi*(diam_in**2))

        self.P_tank_prev = self.P_tank
        self.V_dot_liq_prev = -1e-9 #NOTE: close to zero, but off zero

        self.T_atm = self.T_liq


        ### Debugging code for adiabatic liq, gas and inj liq nodes:
        
        h_ig_gas = analytical_integration_ig_enthalpy(T_REF, self.T_gas)
        u_ig_gas = analytical_integration_ig_int_energy(T_REF, self.T_gas)
        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=self.T_gas, P=P_tank)
        
        h_ig_liq = analytical_integration_ig_enthalpy(T_REF, self.T_liq)
        u_ig_liq = analytical_integration_ig_int_energy(T_REF, self.T_liq)
        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=self.T_liq, P=P_tank)

        self.h_liq_prev = preos_l.H_dep_l/MW + h_ig_liq
        self.h_gas_prev = preos_g.H_dep_g/MW + h_ig_gas

        self.u_liq_prev = preos_l.U_dep_l/MW + u_ig_liq
        self.u_gas_prev = preos_g.U_dep_g/MW + u_ig_gas


        ### below is debug code delete later!!!!

        self.m_inj = 0
        u_liq = preos_l.U_dep_l/MW + u_ig_liq
        u_gas = preos_g.U_dep_g/MW + u_ig_gas
        u_inj = u_liq #NOTE:  this is wrong!!!!

        self.U_liq = u_liq*self.m_liq
        self.U_gas = u_gas*self.m_gas
        self.U_inj = u_inj*self.m_inj

        self.Q_net_liq_cv_out = 0
        self.Q_net_gas_cv_out = 0

        print("starting densities! ", self.rho_liq, self.rho_gas)





    def system_of_liq_odes(self, t, y, P_cc):

        #NOTE: EMPIRICAL FACTOR E:
        E = 2.1e4 #FOR HEAT TRANSFER

        T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas, a, b, c, d, e, f = y  # Unpack state variables

        ### Solve thermo parameters!
        rho_liq, rho_gas, P_tank = solve_thermo_params(T_liq, T_gas, m_liq, m_gas, self.P_tank_prev, self.V_tank, self.volume_err_tolerance)

        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
        T_sat = preos_g.Tsat(P_tank)


        # Mass transfer (1) from injector
        m_dot_inj = spi_model(self.Cd_1, self.A_inj_1, P_tank, P_cc, rho_liq)


        # Heat transfer (2) from saturated surface to gas                       (T_1, T_2, P_tank, rho_2, c, n, tank_diam, fluid)
        # L = tank inner diam , Area of circle x section
        T_film_gas = ((T_sat + T_gas)/2 )
        Q_dot_sat_surf_to_gas = solve_Q_dot_natural_convection_gas(T_sat, T_gas, T_film_gas, P_tank, rho_gas, 0.15, 0.333, self.diam_in, (0.25*np.pi*(self.diam_in**2)), "N2O" ) #relative to gas cv
    
        # Heat transfer (3)  from liq to saturated surface (sat surface assumed to be a liquid with quality 0)
        T_film_liq = ((T_sat + T_liq)/2 )
        Q_dot_liq_to_sat_surf = (E)*solve_Q_dot_natural_convection_liq(T_liq, T_sat, T_film_liq, P_tank, rho_liq, 0.15, 0.333, self.diam_in, (0.25*np.pi*(self.diam_in**2)), "N2O" ) #relative to liq cv
        #NOTE:CORRECTION FACTOR for nitrous oxide heat transfer E = (E) to account for blowing as per [7],[8]

        # Mass transfer (3) by evaporation 
        m_dot_evap = solve_m_dot_evap( T_gas, T_liq, P_tank, Q_dot_liq_to_sat_surf, Q_dot_sat_surf_to_gas)
        #print("m_dot_evap 1: ", m_dot_evap)
        #print("m_dot_evap! ", m_dot_evap, (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas), Q_dot_liq_to_sat_surf, " - ", Q_dot_sat_surf_to_gas)

        # Mass transfer (2) by condensation
        V_gas = m_gas/rho_gas
        V_liq = self.V_tank - V_gas

        m_dot_cond = solve_m_dot_condensed(T_gas, T_liq, P_tank, V_gas, t)

        #preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
        #P_tank = preos_g.Psat(T_gas)


        # Net Mass Transfer of Liquid and Gas CV
        m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj)
        
        #then solve the height of the gas wall
        h_gas_wall = V_gas / (0.25*np.pi*(self.diam_in**2))
        V_gas_wall = 0.25*np.pi*((self.diam_out**2)-(self.diam_in**2))*h_gas_wall
        m_gas_wall = self.rho_wall*V_gas_wall

        h_liq_wall = self.height_tank - h_gas_wall
        V_liq_wall = 0.25*np.pi*((self.diam_out**2)-(self.diam_in**2))*h_liq_wall
        m_liq_wall = self.rho_wall*V_liq_wall

        
        # Heat transfer (4) [natural convection] from liq wall to liq
        Q_dot_liq_wall_to_liq = 0#solve_Q_dot_natural_convection_liq(T_wall_liq, T_liq, T_liq, P_tank, rho_liq, 0.021, 0.4, h_liq_wall, (np.pi*self.diam_in*h_liq_wall), "N2O" ) #relative to liq cv       
        # Heat transfer (5) [natural convection] from gas to gas wall
        Q_dot_gas_wall_to_gas = 0#solve_Q_dot_natural_convection_gas(T_wall_gas, T_gas, T_gas, P_tank, rho_gas, 0.021, 0.4, h_gas_wall, (np.pi*self.diam_in*h_gas_wall), "N2O" ) #relative to gas cv
        

        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
        T_sat = preos_l.Tsat(P_tank)
        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
        preos_sat = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_sat, P=P_tank)

        h_ig_liq = analytical_integration_ig_enthalpy(T_REF, T_liq)
        h_liq = preos_l.H_dep_l/MW + h_ig_liq + H_OFFSET

        h_ig_sat = analytical_integration_ig_enthalpy(T_REF, T_sat)
        h_sat_gas = preos_sat.H_dep_g/MW + h_ig_sat + H_OFFSET
        h_sat_liq = preos_sat.H_dep_l/MW + h_ig_sat + H_OFFSET

        h_ig_gas = analytical_integration_ig_enthalpy(T_REF, T_gas)
        h_gas = preos_g.H_dep_g/MW + h_ig_gas + H_OFFSET

###NOTE: Q_dot eqns: 
        Q_dot_liq = Q_dot_liq_wall_to_liq - Q_dot_liq_to_sat_surf + m_dot_evap*(h_liq - h_sat_liq)
        Q_dot_gas = Q_dot_gas_wall_to_gas + Q_dot_sat_surf_to_gas + m_dot_evap*(h_sat_gas - h_gas)

        #print("Q_dot_liq signs: ", Q_dot_liq, Q_dot_liq_wall_to_liq, - Q_dot_liq_to_sat_surf, m_dot_evap*(h_liq - h_sat_liq) )
        #print("Q_dot_gas signs: ", Q_dot_gas, Q_dot_gas_wall_to_gas,   Q_dot_sat_surf_to_gas, m_dot_evap*(h_sat_gas - h_gas) )
        
        #I believe above is correct but this incase i need to change where terms are applied:
        #Q_dot_liq = Q_dot_liq_wall_to_liq - Q_dot_liq_to_sat_surf + m_dot_evap*( preos_g.Hvap(T_gas)/MW  )  #- m_dot_evap*( preos_g.Hvap(T_gas)/MW  ) I think this is double counting
        #Q_dot_gas = Q_dot_gas_wall_to_gas + Q_dot_sat_surf_to_gas - m_dot_cond*( (-1)*preos_l.Hvap(T_liq)/MW )


        #NOTE: USE ambient properties for air, T_2 will be respective wall temperature (RK var)
        # (6) [natural convection] from atm to liq wall
        Q_dot_atm_to_liq_wall = solve_Q_dot_natural_convection_gas(self.T_atm, T_wall_liq, self.T_atm, self.P_atm, self.rho_atm, 0.59, 0.25, self.height_tank, (np.pi*self.diam_in*h_liq_wall), "Air") #relative to wall_liq cv
        # (7) [natural convection] from atm to gas wall
        Q_dot_atm_to_gas_wall = solve_Q_dot_natural_convection_gas(self.T_atm, T_wall_gas, self.T_atm, self.P_atm, self.rho_atm, 0.59, 0.25, self.height_tank, (np.pi*self.diam_in*h_gas_wall), "Air") #relative to wall_gas cv
        # (8) [conduction] from liq wall to gas wall 
        Q_dot_liq_wall_to_gas_wall = solve_Q_dot_conduction( (T_wall_liq-T_wall_gas), self.height_tank, self.k_w, self.diam_in, self.diam_out) #relative to wall_liq cv



        # Iteratively solve change in CV Volume
        V_dot_liq = self.V_dot_liq_prev #initial guess for dV_dt_liq
        while np.abs(P_dot_error(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, self.h_liq_prev, self.h_gas_prev, self.u_liq_prev, self.u_gas_prev)) > self.P_dot_err_tolerance:
            V_dot_liq = secant((lambda V_dot: P_dot_error(V_dot, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, self.h_liq_prev, self.h_gas_prev, self.u_liq_prev, self.u_gas_prev)), V_dot_liq)
        


        ### Wall nodes:
        
        height_dot = V_dot_liq / (0.25*np.pi*(self.diam_in**2))

        m_dot_liq_wall = self.rho_wall*(0.25*np.pi*height_dot*((self.diam_out**2)-(self.diam_in**2)))  #BUG: this might be a bit unstable w runge kutta steps?
        m_dot_gas_wall = -m_dot_liq_wall

        #NOTE: for T_dot_wall_liq:  IN: (6)      OUT: (4) AND (8)
        T_dot_wall_liq = 0#( Q_dot_atm_to_liq_wall - Q_dot_liq_wall_to_liq - Q_dot_liq_wall_to_gas_wall + m_dot_liq_wall*(0.15)*(T_wall_gas - T_wall_liq) ) / (0.15*m_liq_wall)

        #NOTE: for T_dot_wall_vap:  IN: (7) and (8)      OUT: (5)
        T_dot_wall_gas = 0#( Q_dot_atm_to_gas_wall - Q_dot_gas_wall_to_gas + Q_dot_liq_wall_to_gas_wall + m_dot_gas_wall*(0.15)*( T_wall_liq - T_wall_gas) ) / (0.15*m_gas_wall)


        T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, self.h_liq_prev, self.h_gas_prev, self.u_liq_prev, self.u_gas_prev, True)





        ### Debug Code:
        
        U_dot_liq = solve_U_dot_liq(T_liq, T_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq, Q_dot_liq, self.h_liq_prev, self.h_gas_prev)
        U_dot_gas = solve_U_dot_gas(T_liq, T_gas, P_tank, m_dot_evap, m_dot_cond, (-V_dot_liq), Q_dot_gas, self.h_liq_prev, self.h_gas_prev)

        


        #print("cons energy for adiabatic nodes w liq exit at inj: ", U_dot_liq + U_dot_gas, " = ", m_dot_inj*h_liq )
        #print("cons energy for adiabatic tank w liq exit at inj: ", U_dot_liq + U_dot_gas, " = ", m_dot_inj*h_liq, "expecting 0: ",  (U_dot_liq + U_dot_gas)-(m_dot_inj*h_liq), (P_tank*(V_dot_liq+(-V_dot_liq))) )

        ### solving U_dot_inj:
        
        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
        T_sat = preos_l.Tsat(P_tank)


        h_ig_liq = analytical_integration_ig_enthalpy(T_REF, T_liq)
        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)

        h_liq = h_ig_liq + preos_l.H_dep_l/MW + H_OFFSET

        #print("m_dot_evap 2: ", m_dot_evap)

        return [T_dot_liq, T_dot_gas, m_dot_liq, m_dot_gas, T_dot_wall_liq, T_dot_wall_gas, np.abs(m_dot_inj), U_dot_liq, U_dot_gas, 1, Q_dot_liq, Q_dot_gas]

    def inst(self, P_cc):

        t = 0
        y0 = [self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas, self.m_inj, self.U_liq, self.U_gas, self.U_inj, self.Q_net_liq_cv_out, self.Q_net_gas_cv_out]
        constants = [P_cc, self.h_liq_prev, self.h_gas_prev, self.u_liq_prev, self.u_gas_prev, self.Q_net_liq_cv_out, self.Q_net_gas_cv_out] #constant over 4 rk steps

        #print(" * * * masses: ", self.m_liq, self.m_gas, self.m_inj)

        # NOTE: y is a vector, k1-k4 are derivatives of sol!
        k1 = self.system_of_liq_odes(t, y0, constants)

        y_2 = [y_i + self.TIMESTEP * k1_i / 2 for y_i, k1_i in zip(y0, k1)]
        k2 = self.system_of_liq_odes(t + self.TIMESTEP / 2, y_2, P_cc)

        y_3 = [y_i + self.TIMESTEP * k2_i / 2 for y_i, k2_i in zip(y0, k2)]
        k3 = self.system_of_liq_odes(t + self.TIMESTEP / 2, y_3, P_cc)

        y_4 = [y_i + self.TIMESTEP * k3_i for y_i, k3_i in zip(y0, k3)]
        k4 = self.system_of_liq_odes(t + self.TIMESTEP, y_4, P_cc)

        y =[ y_i + (self.TIMESTEP / 6) * (k1_i + 2*k2_i + 2*k3_i + k4_i) for y_i, k1_i,k2_i,k3_i,k4_i in zip(y0, k1,k2,k3,k4)]

        self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas, self.m_inj, self.U_liq, self.U_gas, self.U_inj, self.Q_net_liq_cv_out, self.Q_net_gas_cv_out = y

        # (4) iteratively solve P_tank to update thermodynamic properties in each node
        #NOTE: this is just to update vals for downstream graphs
        self.rho_liq, self.rho_gas, self.P_tank = solve_thermo_params(self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.P_tank_prev, self.V_tank, self.volume_err_tolerance)

        #update stored vals for RK est and volumes
        self.P_tank_prev = self.P_tank
        self.V_dot_liq_prev = self.V_liq - self.m_liq/self.rho_liq
        self.V_liq = self.m_liq/self.rho_liq
        self.V_gas = self.V_tank - self.V_liq



t = 0
TIMESTEP = 1e-3

P_atm = 1e5 #Pa
T_atm = 286.5 #K
rho_atm = 1.225 #kg/m^3


### Karabeyoglu test case inputs ~ don't know injector used well!

m_nos = 20 #kg
P_tank = 45e5 #Pa
V_tank = 0.0354 #m^3

diam_out = 0.230 #m #NOTE: thesis didn't provide tank geometry, estimated based off of G type nos dimensions (approx equivalent nos mass to Karabeyoglu run tank)
diam_in = 0.450 #m
rho_wall = 2770 #kg/m^3
k_w = 237 #W/(m K)

Cd_1 = 0.425
A_inj_1 = 0.00003 #m^3 NOTE: GUESS
P_cc = 1.03e6 #Pa



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
#def __init__(self, TIMESTEP, T_atm, m_nos, Cd_1, A_inj_1, V_tank, diam_out, diam_in, rho_wall, k_w, P_tank, P_atm, inj_model)
volume_err_tolerance = 1e-8
P_dot_err_tolerance = 10
tank = model(TIMESTEP, m_nos, P_tank, P_cc, P_atm, T_atm, rho_atm, V_tank, diam_out, diam_in, rho_wall, k_w, Cd_1, A_inj_1, volume_err_tolerance, P_dot_err_tolerance)

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

### Debugging liq gas and inj nodes:
init_m_liq = tank.m_liq
init_m_gas = tank.m_gas
init_m_inj = tank.m_inj
init_U_liq = tank.U_liq
init_U_gas = tank.U_gas
init_U_inj = 0

###NOTE:
#is thermo table or thermo lookup faster?

###TODO: try solving single solve different ways!
try:
    start_time = time.time()  # Start timer

    while(t < 500*TIMESTEP): #3000*TIMESTEP
        
        tank.inst(P_cc)
        t+=TIMESTEP 
        LOOKUP_TIME = t

        time_arr.append(t)
        #print("\n next timestep \n")
        
        P_tank_arr.append(tank.P_tank)
        m_tank_arr.append( (tank.m_liq+tank.m_gas+tank.m_inj) )
        m_liq_arr.append(tank.m_liq)
        m_gas_arr.append(tank.m_gas)




        T_liq_arr.append(tank.T_liq)
        T_gas_arr.append(tank.T_gas)

        T_liq_wall_arr.append(tank.T_wall_liq)
        T_gas_wall_arr.append(tank.T_wall_gas)

        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=tank.T_gas, P=tank.P_tank)
        
        T_sat = preos_g.Tsat(tank.P_tank)
        T_sat_arr.append(T_sat)

        P_sat_g = preos_g.Psat(tank.T_gas)
        P_sat_gas_arr.append(P_sat_g)

        P_sat_l = preos_g.Psat(tank.T_liq)
        P_sat_liq_arr.append(P_sat_l)


        ### Debugging Adiabatic liquid gas and inj nodes!!!!
        m_inj_arr.append(tank.m_inj)
        U_liq_arr.append(tank.U_liq)
        U_gas_arr.append(tank.U_gas)


        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=tank.T_liq, P=P_tank)

        u_ig_liq = analytical_integration_ig_int_energy(T_REF, tank.T_liq)
        u_liq = preos_l.U_dep_l/MW + u_ig_liq 

        U_inj_arr.append(u_liq*tank.m_inj)
        #print("u inj sign convention!!! ", tank.u_inj)


        V_liq_arr.append(tank.m_liq/tank.rho_liq)
        V_gas_arr.append(tank.m_gas/tank.rho_gas)

        rho_gas_arr.append(tank.rho_gas)
        rho_liq_arr.append(tank.rho_liq)





        print(f"at t = {t}, final y: ",tank.T_liq, tank.T_gas, tank.m_liq, tank.m_gas, tank.T_wall_liq, tank.T_wall_gas, "\n")

    end_time = time.time()  # End timer
    print(f"\nTotal simulation time: {end_time - start_time:.3f} seconds")

except Exception as e:
    traceback.print_exc()



plt.subplot(1,3,1)
plt.scatter(time_arr,P_tank_arr,label = "tank")
plt.scatter(time_arr,P_sat_liq_arr,label = "P_sat_liq")
plt.scatter(time_arr,P_sat_gas_arr,label = "P_sat_gas")
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Pa)')
plt.title('Pressure vs. Time')
plt.legend()
plt.grid(True)

plt.subplot(1,3,2)
plt.scatter(time_arr,m_liq_arr, label = "liquid")
plt.scatter(time_arr,m_gas_arr, label = "gas")
plt.scatter(time_arr,m_tank_arr, label = "total mass")
plt.xlabel('Time (s)')
plt.ylabel('Mass (kg)')
plt.title('Mass vs. Time')
plt.legend()
plt.grid(True)

plt.subplot(1,3,3)
plt.scatter(time_arr,T_liq_arr, label = "liquid")
plt.scatter(time_arr,T_gas_arr, label = "gas")
plt.scatter(time_arr,T_sat_arr, label = "T_sat")
plt.scatter(time_arr,T_liq_wall_arr, label = "WALL liquid")
plt.scatter(time_arr,T_gas_wall_arr, label = "WALL gas")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.title('Temperature vs. Time')
plt.legend()
plt.grid(True)
plt.show()


plt.subplot(1,2,1)
plt.scatter(time_arr,rho_liq_arr, label = "liquid")
plt.scatter(time_arr,rho_gas_arr, label = "gas")
plt.xlabel('Time (s)')
plt.ylabel('Density (kg/m^3)')
plt.title('Density vs. Time')
plt.legend()
plt.grid(True)


plt.subplot(1,2,2)
plt.scatter(time_arr,V_liq_arr, label = "liquid")
plt.scatter(time_arr,V_gas_arr, label = "gas")
plt.xlabel('Time (s)')
plt.ylabel('Volume (m^3)')
plt.title('Volume vs. Time')
plt.legend()
plt.grid(True)
plt.show()

### Debugging adiabatic liq gas and inj nodes
plt.subplot(1,2,1)
plt.scatter(time_arr,m_liq_arr, label = "liquid")
plt.scatter(time_arr,m_gas_arr, label = "gas")
plt.scatter(time_arr,m_inj_arr, label = "inj")
plt.xlabel('Time (s)')
plt.ylabel('Mass (kg)')
plt.title('Mass vs. Time')
plt.legend()
plt.grid(True)

plt.subplot(1,2,2)
plt.scatter(time_arr, U_liq_arr, label = "liquid")
plt.scatter(time_arr, U_gas_arr, label = "gas")
plt.scatter(time_arr, U_inj_arr, label = "inj")
plt.xlabel('Time (s)')
plt.ylabel('Energy (J)')
plt.title('Energy vs. Time')
plt.legend()
plt.grid(True)


print(f"\n\n\n\nInitial Total Mass: {init_m_liq + init_m_gas + np.abs(init_m_inj)} (kg)")

percent_diff_m = ((tank.m_liq + tank.m_gas + tank.m_inj)-(init_m_liq + init_m_gas + init_m_inj)) /((init_m_liq + init_m_gas + init_m_inj))*100


print("\nMass of gas before/after: ", init_m_gas, tank.m_gas, "\nMass of liq before/after: ",init_m_liq, tank.m_liq)

plt.show()