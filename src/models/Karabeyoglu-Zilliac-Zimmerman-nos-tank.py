from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical
import matplotlib.pyplot as plt

# Global Constants:
R_U = 8.31446 #J/(mol K)
T_REF = 298.15 #K
P_REF = 101325 #Pa

n2o = Chemical('nitrous oxide', T=T_REF)

MW = (n2o.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o.omega - 0.26992*n2o.omega**2

TANK_DIAM = 0.0254*5.5 #m
CS_AREA = 0.25*np.pi*(TANK_DIAM**2) #m^2
g = 9.81 #m/s^2

#this is failing because we are inputting zero here. (why x2-x1=0)
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
        kk = kk + 1
    x = x2
    return x

#NOTE: SIGN CONVENTION: Q_dot dir: T_1 --> T_2 (_2 is the fluid transferring to!)
def solve_Q_dot_natural_convection_liq(T_1, T_2, P_tank, rho_2, c, n, tank_diam, fluid): #BUG: potential mistake, solving _1 properties with _2 inputs, double check this is likely a mistake

    n2o = Chemical('N2O', T=T_2, P=P_tank) 
    k_2 = n2o.kl
    visc_2 = n2o.mul

    preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_2, P=P_tank)
    Cp_2 = preos_l.Cp_dep_l + n2o.Cpg 

    dV_dT_P = preos_l.dV_dT_l

    beta = dV_dT_P*rho_2 

    Gr = ((tank_diam**3)*(rho_2**2)*g*beta*np.abs(T_2 - T_1) ) / (visc_2**2)
    Pr = (Cp_2*visc_2)/ k_2

    X = Gr*Pr
    h = c * (k_2/tank_diam) * X**n

    Q_dot = h*(0.25*np.pi*(tank_diam**2))*(T_1-T_2)

    return Q_dot #NOTE: Q_dot + going into (2)

#NOTE: SIGN CONVENTION: Q_dot dir: T_1 --> T_2 (_2 is the fluid sign relative to!)
def solve_Q_dot_natural_convection_gas(T_1, T_2, P_tank, rho_2, c, n, tank_diam, fluid): #BUG: potential mistake, solving _1 properties with _2 inputs, double check this is likely a mistake

    n2o = Chemical('N2O', T=T_2, P=P_tank)  #TODO: units here!!!
    k_2 = n2o.kg
    visc_2 = n2o.mug

    preos_g = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_2, P=P_tank)
    Cp_2 = preos_g.Cp_dep_g + n2o.Cpg

    dV_dT_P = preos_g.dV_dT_g

    beta = dV_dT_P*rho_2 

    Gr = ((tank_diam**3)*(rho_2**2)*g*beta*np.abs(T_2 - T_1) ) / (visc_2**2)
    Pr = (Cp_2*visc_2)/ k_2

    X = Gr*Pr
    h = c * (k_2/tank_diam) * X**n

    Q_dot = h*(0.25*np.pi*(tank_diam**2))*(T_1-T_2)

    return Q_dot


def solve_Q_dot_conduction(delta_T, h_tank, k_w, diam_in, diam_out):

    L_w_cond = 0.5*h_tank
    Q_dot_conduction = k_w *(delta_T)*(0.25*np.pi*((diam_out**2)-(diam_in**2)))/L_w_cond
    
    return Q_dot_conduction

#TODO: update with other injector model once we get this thing up
def spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit):
    m_dot_spi = Cd_hem_spi_dyer * A_inj_ox * np.sqrt( 2 * rho_tank_exit * (P_1 - P_2)  )
    return m_dot_spi

def solve_m_dot_condensed(T_gas, P, V_gas, t):

    preos = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_gas, P=P_tank)
    P_sat = preos.Psat(T_gas)

    if P > P_sat:
        return ((P-P_sat)*V_gas*MW)/( (R_U/MW)*T_gas*t)
    else:
        return 0


def solve_U_dot_liq(T_liq, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq , Q_dot_net):
    n2o_ig = Chemical('N2O', T=T_liq) 
    preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_liq, P=P_tank)
    h_liq = preos_l.H_dep_l + (n2o_ig.Cpg*T_liq - n2o.Cpg*T_REF) + R_U*T_REF/MW
    latent_heat_evap = preos_l.Hvap(T_liq)
    U_dot = m_dot_inj*h_liq + m_dot_evap*latent_heat_evap + m_dot_cond*(-1*latent_heat_evap) - P_tank*V_dot_liq + Q_dot_net #NOTE: (-1)* to convert latent heat evap to condensation
    return U_dot

def solve_U_dot_gas(T_liq, P_tank, m_dot_evap, m_dot_cond, V_dot_liq , Q_dot_net):
    preos_g = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_liq, P=P_tank)
    latent_heat_evap = preos_g.Hvap(T_liq)
    U_dot = m_dot_evap*latent_heat_evap + m_dot_cond*(-1*latent_heat_evap) - P_tank*V_dot_liq + Q_dot_net #NOTE: (-1)* to convert latent heat evap to condensation
    return U_dot


#BUG: NOW PRESSURE IS NEGATIVE HERE?
def V_tank_error(P_guess, T_liq, T_gas, m_liq, m_gas, V_tank):

    preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_liq, P=P_guess)
    rho_liq = preos_l.rho_l*MW
    
    preos_g = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_gas, P=P_guess)
    rho_gas = preos_g.rho_g*MW

    V_tank_est = (m_liq/rho_liq) + (m_gas/rho_gas)

    #print("   V_tank_err: ", V_tank_est - V_tank )

    return V_tank_est - V_tank 

def solve_thermo_params(T_liq, T_gas, m_liq, m_gas, P_tank_prev, V_tank, volume_err_tolerance):

    P_tank = P_tank_prev #initial guess for pressure

    while np.abs(V_tank_error(P_tank, T_liq, T_gas, m_liq, m_gas, V_tank) ) > volume_err_tolerance:
        P_tank = secant((lambda P: V_tank_error(P, T_liq, T_gas, m_liq, m_gas, V_tank)), P_tank)

    #print("solved P_tank, ", P_tank, T_liq, T_gas)

    preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_liq, P=P_tank)
    rho_liq = preos_l.rho_l*MW
    
    preos_g = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_gas, P=P_tank)
    rho_gas = preos_g.rho_g*MW

    return rho_liq, rho_gas, P_tank

def single_solve_T_dot_liq_gas(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas):

    m_dot_liq = -m_dot_evap + m_dot_cond - m_dot_inj
    m_dot_gas = m_dot_evap - m_dot_cond #convert sign convention from liq cv to gas cv

    V_dot_gas = -V_dot_liq

    V_liq = m_liq/rho_liq
    V_gas = m_gas/rho_gas

    rho_dot_liq = (1 / V_liq) * m_dot_liq - (m_liq / V_liq**2) * V_dot_liq
    rho_dot_gas = (1 / V_gas) * m_dot_gas - (m_gas / V_gas**2) * V_dot_gas  


    #(T_liq, P_tank, m_dot_evap, m_dot_cond, V_dot_liq , Q_dot_net)
    U_dot_liq = solve_U_dot_liq(T_liq, P_tank, (-1)*m_dot_inj, (-1)*m_dot_evap, m_dot_cond, V_dot_liq , Q_dot_liq)
    U_dot_gas = solve_U_dot_gas(T_gas, P_tank, m_dot_evap, (-1)*m_dot_cond, V_dot_gas, Q_dot_gas)

    preos_g = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_gas, P=P_tank)
    partial_du_dv_const_T_gas = T_gas*preos_g.dP_dT_g - P_tank
    #print("\n partial: ", partial_du_drho_const_T_gas, "\n")
    n2o_ig = Chemical('N2O', T=T_gas) 
    Cv_gas = n2o_ig.Cvg + preos_g.Cv_dep_g
    u_gas = preos_g.U_dep_g/MW + (n2o_ig.Cvg*T_gas - n2o.Cvg*T_REF)

    preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_liq, P=P_tank)
    partial_du_dv_const_T_liq = T_liq*preos_l.dP_dT_l - P_tank 
    Cv_liq = n2o_ig.Cvg + preos_l.Cv_dep_l #NOTE: this might be wrong but i think with departure it is supposed to be referencing ideal gas property
    u_liq = preos_l.U_dep_l/MW + (n2o_ig.Cvg*T_gas - n2o.Cvg*T_REF) #NOTE: same warning as ^ but also for this line

    #print( Cv_gas, m_gas, U_dot_gas, u_gas, m_dot_gas, partial_du_dv_const_T_gas, rho_dot_gas)
   
    T_dot_liq = (1/Cv_liq)*( (1/m_liq) * (U_dot_liq - u_liq*m_dot_liq) - (partial_du_dv_const_T_liq* V_dot_liq) )
    T_dot_gas = (1/Cv_gas)*( (1/m_gas) * (U_dot_gas - u_gas*m_dot_gas) - (partial_du_dv_const_T_gas* V_dot_gas) )
    
    return T_dot_liq, T_dot_gas #TODO: TRACE THROUGH T_DOT SOLVING

def P_dot_error(V_dot_guess, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas):    

    V_liq = m_liq/rho_liq
    V_gas = m_gas/rho_gas

    V_dot_gas = -V_dot_guess

    T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_guess, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas)

    #print("  P_dot_error: ",T_dot_liq, T_dot_gas)

    preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_liq, P=P_tank)
    partial_dP_dT_const_v_liq = preos_l.dP_dT_l
    partial_dP_dv_const_T_liq = preos_l.dP_dV_l

    P_dot_liq = partial_dP_dT_const_v_liq*T_dot_liq + partial_dP_dv_const_T_liq*V_dot_guess

    preos_g = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_gas, P=P_tank)
    partial_dP_dT_const_v_gas = preos_g.dP_dT_g
    partial_dP_dv_const_T_gas = preos_g.dP_dV_g

    P_dot_gas = partial_dP_dT_const_v_gas*T_dot_gas + partial_dP_dv_const_T_gas*V_dot_gas

    return P_dot_liq - P_dot_gas

def iteratively_solve_T_dot_liq_and_gas(m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, V_dot_liq_prev, P_dot_err_tolerance):

    V_dot_liq = V_dot_liq_prev #initial guess for dV_dt_liq

    while np.abs(P_dot_error(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas) ) > P_dot_err_tolerance:
        V_dot_liq = secant((lambda V_dot: P_dot_error(V_dot, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas )), V_dot_liq)

    T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas)

    #print("  iter solve T_dots: ",T_dot_liq, T_dot_gas)
    return T_dot_liq, T_dot_gas


def system_of_liq_odes(t, y, constants):

    T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas = y  # Unpack state variables

    T_atm, P_atm, rho_atm, height_tank, k_w, diam_in, diam_out, Cd_1, A_inj_1, P_cc, V_tank, rho_wall, P_tank_prev, V_dot_liq_prev, volume_err_tolerance, P_dot_err_tolerance = constants

    rho_liq, rho_gas, P_tank = solve_thermo_params(T_liq, T_gas, m_liq, m_gas, P_tank_prev, V_tank, volume_err_tolerance)

    preos = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_gas, P=P_tank)
    T_sat = preos.Tsat(P_tank)

    #print(rho_liq, rho_gas, P_tank)

    # (2) from saturated surface to gas                       (T_1, T_2, P_tank, rho_2, c, n, tank_diam, fluid)
    Q_dot_sat_surf_to_gas = solve_Q_dot_natural_convection_gas(T_sat, T_gas, P_tank, rho_gas, 0.15, 0.333, diam_in, "N2O") #relative to gas cv
    # (3)  from liq to saturated surface (sat surface assumed to be a liquid with quality 0)
    Q_dot_liq_to_sat_surf = (-1)* solve_Q_dot_natural_convection_liq(T_sat, T_liq, P_tank, rho_liq, 0.15, 0.333, diam_in, "N2O") #relative to liq cv
    # (4) [natural convection] from liq wall to liq
    Q_dot_liq_wall_to_liq = solve_Q_dot_natural_convection_liq(T_wall_liq, T_liq, P_tank, rho_liq, 0.59, 0.25, diam_in, "N2O") #relative to liq cv       
    # (5) [natural convection] from gas wall to gas
    Q_dot_gas_wall_to_gas = solve_Q_dot_natural_convection_gas(T_wall_gas, T_gas, P_tank, rho_gas, 0.59, 0.25, diam_in, "N2O") #relative to gas cv

    #NOTE: USE ambient properties for air, T_2 will be respective wall temperature (RK var)
    # (6) [natural convection] from atm to liq wall
    Q_dot_atm_to_liq_wall = (-1)* solve_Q_dot_natural_convection_gas(T_wall_liq, T_atm, P_atm, rho_atm, 0.59, 0.25, diam_out, "air") #relative to wall_liq
    # (7) [natural convection] from atm to gas wall
    Q_dot_atm_to_gas_wall = (-1)* solve_Q_dot_natural_convection_gas(T_wall_gas, T_atm, P_atm, rho_atm, 0.59, 0.25, diam_out, "air") #relative to wall_gas
    # (8) [conduction] from liq wall to gas wall
    Q_dot_liq_wall_to_gas_wall = solve_Q_dot_conduction( (T_wall_gas-T_wall_liq), height_tank, k_w, diam_in, diam_out) #relative to liquid cv

    #print( Q_dot_atm_to_liq_wall, Q_dot_atm_to_gas_wall, Q_dot_liq_wall_to_gas_wall) #NOTE: in thermal eq, will need to check next rk step once solve this issue

    #(1) mass transfer from injector already solved
    m_dot_inj = spi_model(Cd_1, A_inj_1, P_tank, P_cc, rho_liq)
    #(2) mass transfer by condensation
    V_gas = m_gas*rho_gas
    m_dot_cond = solve_m_dot_condensed(T_gas, P_tank, V_gas, t)
    #(3) mass transfer by evaporation

    preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_liq, P=P_tank)
    latent_heat_evap = preos_l.Hvap(T_liq)

    m_dot_evap = (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) / latent_heat_evap #NOTE: P sure its T_liq but might be wrong

    m_dot_liq = -m_dot_evap + m_dot_cond - m_dot_inj
    m_dot_gas = m_dot_evap - m_dot_cond #convert sign convention from liq cv to gas cv

    ####solve wall properties!!! then heat transfer terms

    #TODO: solve m_dot_wall

    Q_dot_liq = -Q_dot_liq_wall_to_liq + Q_dot_liq_to_sat_surf
    Q_dot_gas = Q_dot_gas_wall_to_gas*(-1)
    #print("\nchecking sign convention Q_dot: ", Q_dot_liq, Q_dot_gas)

    V_dot_liq = V_dot_liq_prev #initial guess for dV_dt_liq

    while np.abs(P_dot_error(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas) ) > P_dot_err_tolerance:
        V_dot_liq = secant((lambda V_dot: P_dot_error(V_dot, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas )), V_dot_liq)

    delta_height = V_dot_liq / (0.25*np.pi*(diam_in**2))

    m_dot_liq_wall = 0
    if t != 0: m_dot_liq_wall = rho_wall*(0.25*np.pi*delta_height*((diam_out**2)-(diam_in**2)))/t
    m_dot_gas_wall = - m_dot_liq_wall

    m_gas_wall = V_gas / (0.25*np.pi*(diam_in**2))
    m_liq_wall = (V_tank-V_gas) / (0.25*np.pi*(diam_in**2))

    
    #NOTE: for T_dot_wall_liq:  IN: (6)      OUT: (4) AND (8)
    T_dot_wall_liq = ( Q_dot_atm_to_liq_wall + Q_dot_liq_wall_to_liq - Q_dot_liq_wall_to_gas_wall +  m_dot_liq_wall*(0.15)*( T_wall_gas - T_wall_liq) ) / (0.15*m_liq_wall)

    #NOTE: for T_dot_wall_vap:  IN: (7) and (8)      OUT: (5)
    T_dot_wall_gas = ( -Q_dot_atm_to_gas_wall + Q_dot_gas_wall_to_gas + Q_dot_liq_wall_to_gas_wall +  m_dot_gas_wall*(0.15)*(T_wall_liq - T_wall_gas) ) / (0.15*m_gas_wall)



    ### iteratively solving dV/dt with pressure constraint, dT/dt --> f( dU/dt and d_rho/dt --> f( dV/dt) )
    T_dot_liq, T_dot_gas = iteratively_solve_T_dot_liq_and_gas(m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, V_dot_liq_prev, P_dot_err_tolerance)

    return [T_dot_liq, T_dot_gas, m_dot_liq, m_dot_gas, T_dot_wall_liq, T_dot_wall_gas]




class model():

    def __init__(self, timestep, m_nos, P_tank, P_cc, P_atm, T_atm, rho_atm, V_tank, diam_out, diam_in, rho_wall, k_w, Cd_1, A_inj_1, volume_err_tolerance, P_dot_err_tolerance):
        
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

        #NOTE: Setup Thermo Properties - Assuming Tank Starts NEAR Sat Conditions

        self.V_tank = V_tank
        self.m_nos = m_nos

        preos = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=self.T_atm, P=self.P_tank)
        rho_sat_gas = (1/preos.V_g_sat(self.T_atm))*MW
        rho_sat_liq = (1/preos.V_l_sat(self.T_atm))*MW
        
        rho_sat_tank = self.m_nos/V_tank

        x_tank = ( (1/rho_sat_tank)-(1/rho_sat_liq) ) / ( (1/rho_sat_gas)-(1/rho_sat_liq) )

        #gas cv
        self.m_gas = x_tank*self.m_nos

        self.T_gas = preos.Tsat(self.P_tank) + 0.05 #Perturb to start close to equillibrium but as a gas
        self.T_wall_gas = self.T_gas

        preos = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=self.T_gas, P=self.P_tank)
        self.rho_gas = preos.rho_g*MW
        v_gas = 1/self.rho_gas
        self.V_gas = v_gas*self.m_gas

        #liquid cv
        self.m_liq = self.m_nos-self.m_gas

        preos = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=self.T_gas, P=self.P_tank)
        self.T_liq = preos.Tsat(self.P_tank) - 0.05 #Perturb to start close to equillibrium
        self.T_wall_liq = self.T_liq

        #solve rho_liq!
        preos = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=self.T_liq, P=self.P_tank)
        self.rho_liq = preos.rho_l*MW
        v_liq = 1/self.rho_liq
        self.V_liq = v_liq*self.m_liq

        #self.rho_exit = 1/v_liq

        self.V_tank = self.V_liq+self.V_gas # "what are you going to do if the aluminum is too small? water it? give it sunlight? let it grow?"
        self.height_tank = self.V_tank/(0.25*np.pi*(diam_in**2))
        self.P_tank_prev = self.P_tank
        self.V_dot_liq_prev = 1e-6 #NOTE: close to zero, but off zero


    def inst(self, P_cc):
        
        #TODO: check for liquid phase!!!!

        # (1) enter Explicit Runge Kutta 6 ode system
        #T_atm, P_atm, rho_atm, height_tank, k_w, diam_in, diam_out, Cd_1, A_inj_1, P_cc, V_tank, rho_wall, P_tank_prev, V_dot_liq_prev, volume_err_tolerance, P_dot_err_tolerance
        constants = [self.T_atm, self.P_atm, self.rho_atm, self.height_tank, self.k_w, self.diam_in, self.diam_out, self.Cd_1, self.A_inj_1, self.P_cc, self.V_tank, self.rho_wall, self.P_tank_prev, self.V_dot_liq_prev, self.volume_err_tolerance, self.P_dot_err_tolerance]
        y0 =  [self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas ]

        #solve with scipy
        solution = solve_ivp(system_of_liq_odes, [0, self.timestep], y0, args=(constants,), method='RK45') #, rtol=1e-12, atol=1e-12) #NOTE: THIS IS PROBABLY TOO SLOW, I SET THIS HIGH TOLERANCE FOR THE LAST TRY
            

        # Extract the final state from the solution
        final_state = solution.y[:, -1]  # Get the last column, which represents the state at the final time step

        # Unpack the final state into your variables
        self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas = final_state

        # (4) iteratively solve P_tank to update thermodynamic properties in each node
        self.rho_liq, self.rho_gas, self.P_tank = solve_thermo_params(self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.P_tank_prev, self.V_tank, self.volume_err_tolerance)

        #update stored vals for RK est and volumes
        self.P_tank_prev = self.P_tank
        self.V_dot_liq_prev = self.V_liq - self.m_liq/self.rho_liq
        self.V_liq = self.m_liq/self.rho_liq
        self.V_gas = self.V_tank - self.V_liq


t = 0
TIMESTEP = 1e-4

P_atm = 1e5 #Pa
T_atm = 273.15 + 15 #K
rho_atm = 1.225 #kg/m^3

m_nos = 20
P_tank = 45e5
V_tank = 0.0354

diam_out = 0.230 #m #NOTE: thesis didn't provide tank geometry, estimated based off of G type nos dimensions (approx equivalent mass to Karabeyoglu run tank)
diam_in = 0.215 #m
rho_wall = 2770 #kg/m^3
k_w = 237 #W/(m K)

Cd_1 = 0.425
A_inj_1 = 0.00003 #m^3 NOTE: GUESS
P_cc = 1.03e6

inj_model = None #TODO: implement

#def __init__(self, TIMESTEP, T_atm, m_nos, Cd_1, A_inj_1, V_tank, diam_out, diam_in, rho_wall, k_w, P_tank, P_atm, inj_model)
volume_err_tolerance = 1e-7
P_dot_err_tolerance = 1e2
tank = model(TIMESTEP, m_nos, P_tank, P_cc, P_atm, T_atm, rho_atm, V_tank, diam_out, diam_in, rho_wall, k_w, Cd_1, A_inj_1, volume_err_tolerance, P_dot_err_tolerance)

v_liq_arr = []
P_tank_arr = []
v_vap_arr = []
time_arr = []
#try:


TIMESTEP = 1e-4


while(t<1000*TIMESTEP):
    tank.inst(P_cc)
    t+=TIMESTEP 

    time_arr.append(t)
    v_liq_arr.append(1/tank.rho_liq)
    P_tank_arr.append(tank.P_tank)
    v_vap_arr.append(1/tank.rho_gas)


print("done")




plt.plot(time_arr,P_tank_arr)
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Pa)')
plt.title('Pressure vs. Time')
plt.grid(True)
plt.legend()
plt.show()