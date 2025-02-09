from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical
import matplotlib.pyplot as plt
import traceback

# Global Constants:
R_U = 8.31446 #J/(mol K)
T_REF = 298.15 #K

n2o_g = Chemical('nitrous oxide', T=T_REF)
PC = n2o_g.Pc
TC = n2o_g.Tc
OMEGA = n2o_g.omega

MW = (n2o_g.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o_g.omega - 0.26992*n2o_g.omega**2

TANK_DIAM = 0.0254*5.5 #m
CS_AREA = 0.25*np.pi*(TANK_DIAM**2) #m^2
g = 9.81 #m/s^2

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

#NOTE: SIGN CONVENTION: Q_dot dir: (+) T_1 --> T_2 (f is fluid)
def solve_Q_dot_natural_convection_liq(T_1, T_2, T_f, P_tank, rho_f, c, n, tank_diam, fluid):

    n2o = Chemical('N2O', T=T_f, P=P_tank) 
    k_f = n2o.kl #W/m/K
    visc_f = n2o.mul #Pa s

    preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_f, P=P_tank)
    Cp_2 = (preos_l.Cp_dep_l/MW + n2o.Cpg) #J/K #NOTE: should be fixed

    dV_dT_P = preos_l.dV_dT_l #BUG! lib documentation "volume"

    beta = dV_dT_P*rho_f 

    Gr = ((tank_diam**3)*(rho_f**2)*g*beta*np.abs(T_2 - T_1) ) / (visc_f**2)
    Pr = (Cp_2*visc_f)/ k_f

    X = Gr*Pr
    h = c * (k_f/tank_diam) * X**n

    Q_dot = h*(0.25*np.pi*(tank_diam**2))*(T_1-T_2)

    return Q_dot #NOTE: Q_dot + going into (2)

#NOTE: SIGN CONVENTION: Q_dot dir: (+) T_1 --> T_2 (f is fluid)
def solve_Q_dot_natural_convection_gas(T_1, T_2, T_f, P_tank, rho_2, c, n, tank_diam, fluid): #BUG: potential mistake, solving _1 properties with _2 inputs, double check this is likely a mistake

    n2o = Chemical('N2O', T=T_f, P=P_tank)  #TODO: units here!!!
    k_f = n2o.kg
    visc_f = n2o.mug

    preos_g = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_f, P=P_tank)
    Cp_2 = (preos_g.Cp_dep_g/MW + n2o.Cpg)

    dV_dT_P = preos_g.dV_dT_g

    beta = dV_dT_P*rho_2 

    Gr = ((tank_diam**3)*(rho_2**2)*g*beta*np.abs(T_2 - T_1) ) / (visc_f**2)
    Pr = (Cp_2*visc_f)/ k_f

    X = Gr*Pr
    h = c * (k_f/tank_diam) * X**n

    Q_dot = h*(0.25*np.pi*(tank_diam**2))*(T_1-T_2)

    return Q_dot


def solve_Q_dot_conduction(delta_T, h_tank, k_w, diam_in, diam_out):

    L_w_cond = 0.3*h_tank
    Q_dot_conduction = k_w *(delta_T)*(0.25*np.pi*((diam_out**2)-(diam_in**2)))/L_w_cond
    
    return Q_dot_conduction

#TODO: update with other injector model once we get this thing up
def spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit):
    m_dot_spi = (-1)*Cd_hem_spi_dyer * A_inj_ox * np.sqrt( 2 * rho_tank_exit * (P_1 - P_2)  )
    return m_dot_spi

def solve_m_dot_condensed(T_gas, P_tank, V_gas, t):
    m_dot_cond = 0

    
    preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
    P_sat = preos_g.Psat(T_gas)

    if (P_tank > P_sat):
        m_dot_cond = ((P_tank-P_sat)*V_gas*MW)/( preos_g.Z_g*(R_U/MW)*T_gas*(t) )  
    
    return m_dot_cond


def solve_U_dot_liq(T_liq, T_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq , Q_dot_net):
    n2o_ig = Chemical('N2O', T=T_liq) 
    preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    h_liq = preos_l.H_dep_l/MW + n2o_ig.H #NOTE: should be fixed now
    latent_heat_evap_l = preos_l.Hvap(T_liq)/MW 

    preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
    latent_heat_cond_g = (-1)*preos_g.Hvap(T_gas)/MW 

    U_dot_liq = -m_dot_inj*h_liq - m_dot_evap*latent_heat_evap_l + m_dot_cond*latent_heat_cond_g - P_tank*V_dot_liq + Q_dot_net

    #print("U_dot_liq: ", U_dot_liq , -m_dot_inj*h_liq ,- m_dot_evap*latent_heat_evap_l ,+ m_dot_cond*(latent_heat_cond_g) , Q_dot_net )

    
    return U_dot_liq

def solve_U_dot_gas(T_liq, T_gas, P_tank, m_dot_evap, m_dot_cond, V_dot_gas, Q_dot_net):
    preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    latent_heat_evap_l = preos_l.Hvap(T_liq)/MW 

    preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
    latent_heat_cond_g = (-1)*preos_g.Hvap(T_gas)/MW 

    U_dot_gas = m_dot_evap*latent_heat_evap_l - m_dot_cond*latent_heat_cond_g - P_tank*V_dot_gas + Q_dot_net 
    
    #print("U_dot_gas: ", U_dot_gas, m_dot_evap*latent_heat_evap_l, - m_dot_cond*(latent_heat_cond_g) , Q_dot_net)
    #print("m_dot_evap! ", m_dot_evap)

    return U_dot_gas

def solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj):
    m_dot_liq = -m_dot_evap + m_dot_cond + m_dot_inj
    m_dot_gas = m_dot_evap - m_dot_cond #convert sign convention from liq cv to gas cv

    return m_dot_liq, m_dot_gas

def V_tank_error(P_guess, T_liq, T_gas, m_liq, m_gas, V_tank):

    #print("this should be a liquid! ", T_liq, P_guess)
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

def single_solve_T_dot_liq_gas(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas):

    m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj)

    V_dot_gas = -V_dot_liq

    d_rho_dt_liq = (1/V_liq)*m_dot_liq -(m_liq/(V_liq**2))*V_dot_liq
    d_rho_dt_gas = (1/V_gas)*m_dot_gas -(m_gas/(V_gas**2))*V_dot_gas

    U_dot_liq = solve_U_dot_liq(T_liq, T_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq, Q_dot_liq)
    U_dot_gas = solve_U_dot_gas(T_liq, T_gas, P_tank, m_dot_evap, m_dot_cond, V_dot_gas, Q_dot_gas) #NOTE: convert evap and cond terms to gas cv


    preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
    partial_du_d_rho_const_T_gas = ((-1)/(rho_gas**2))*(T_gas*preos_g.dP_dT_g - P_tank) #NOTE: deriv wrt Vm -> deriv wrt rho
    n2o_ig_g = Chemical('N2O', T=T_gas) 
    
    Cv_gas = preos_g.Cv_dep_g/MW + n2o_ig_g.Cvg 
    u_gas = preos_g.U_dep_g/MW + (n2o_ig_g.H - (R_U*T_gas/MW) )#NOTE: THIS WAS T_LIQ BUT I THINK IT SHOULD BE T_GAS - no visual difference?

    preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    partial_du_d_rho_const_T_liq = ((-1)/(rho_liq**2))*(T_liq*preos_l.dP_dT_l - P_tank) ##NOTE: deriv wrt Vm -> deriv wrt rho
    n2o_ig_l = Chemical('N2O', T=T_liq) 

    Cv_liq = preos_l.Cv_dep_l/MW + n2o_ig_l.Cvg
    u_liq = preos_l.U_dep_l/MW + (n2o_ig_l.H - (R_U*T_liq/MW) )#NOTE: same warning as ^ but also for this line
   

    T_dot_liq = (1/Cv_liq)*( (1/m_liq) * (U_dot_liq - u_liq*m_dot_liq) - (partial_du_d_rho_const_T_liq* d_rho_dt_liq) )
    T_dot_gas = (1/Cv_gas)*( (1/m_gas) * (U_dot_gas - u_gas*m_dot_gas) - (partial_du_d_rho_const_T_gas* d_rho_dt_gas) )

    print("T_dot_liq! ", T_dot_liq, (1/m_liq) * (U_dot_liq - u_liq*m_dot_liq), -(partial_du_d_rho_const_T_liq* d_rho_dt_liq) )
    print("T_dot_gas! ", T_dot_gas, (1/m_gas) * (U_dot_gas - u_gas*m_dot_gas), -(partial_du_d_rho_const_T_gas* d_rho_dt_gas) )
    #print("T_dot_gas! ", (U_dot_gas - u_gas*m_dot_gas),U_dot_gas, - u_gas*m_dot_gas, u_gas, m_dot_gas)


    return T_dot_liq, T_dot_gas 

def P_dot_error(V_dot_guess, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas):    

    V_dot_gas = -V_dot_guess #guessing for liquid

    m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj) #or put this outside and pass in?

    d_rho_dt_liq = (1/V_liq)*m_dot_liq - (m_liq/(V_liq**2))*V_dot_guess
    d_rho_dt_gas = (1/V_gas)*m_dot_gas - (m_gas/(V_gas**2))*V_dot_gas


    T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_guess, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas)

    preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    partial_dP_dT_const_rho_liq = preos_l.dP_dT_l #NOTE: deriv wrt Vm = deriv wrt rho

    partial_dP_drho_const_T_liq = (preos_l.dP_dV_l)*(-MW/(rho_liq**2)) #converting deriv wrt Vm --> wrt rho

    P_dot_liq = partial_dP_dT_const_rho_liq*T_dot_liq + partial_dP_drho_const_T_liq*d_rho_dt_liq

    preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
    partial_dP_dT_const_rho_gas = preos_g.dP_dT_g #NOTE: deriv wrt Vm = deriv wrt rho

    partial_dP_drho_const_T_gas = (preos_g.dP_dV_g)*(-MW/(rho_gas**2)) #converting deriv wrt Vm --> wrt rho

    P_dot_gas = partial_dP_dT_const_rho_gas*T_dot_gas + partial_dP_drho_const_T_gas*d_rho_dt_gas


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
        self.V_dot_liq_prev = 1e-9 #NOTE: close to zero, but off zero



    def system_of_liq_odes(self, t, y, P_cc):

        T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas = y  # Unpack state variables

        

        rho_liq, rho_gas, P_tank = solve_thermo_params(T_liq, T_gas, m_liq, m_gas, self.P_tank_prev, self.V_tank, self.volume_err_tolerance)

        preos = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
        T_sat = preos.Tsat(P_tank)


        # (2) from saturated surface to gas                       (T_1, T_2, P_tank, rho_2, c, n, tank_diam, fluid)
        Q_dot_sat_surf_to_gas = solve_Q_dot_natural_convection_gas(T_sat, T_gas, T_gas, P_tank, rho_gas, 0.15, 0.333, self.diam_in, "N2O") #relative to gas cv
        # (3)  from liq to saturated surface (sat surface assumed to be a liquid with quality 0)
        Q_dot_liq_to_sat_surf = solve_Q_dot_natural_convection_liq(T_liq, T_sat, T_liq, P_tank, rho_liq, 0.15, 0.333, self.diam_in, "N2O") #relative to liq cv
        #print("Q_dot_liq_to_sat_surf: ", Q_dot_liq_to_sat_surf , T_liq, T_sat, T_liq-T_sat, P_tank, rho_liq)
        
        
        # (4) [natural convection] from liq wall to liq
        Q_dot_liq_wall_to_liq = solve_Q_dot_natural_convection_liq(T_wall_liq, T_liq, T_liq, P_tank, rho_liq, 0.021, 0.4, self.diam_in, "N2O") #relative to liq cv       
        # (5) [natural convection] from gas to gas wall
        Q_dot_gas_wall_to_gas = solve_Q_dot_natural_convection_gas(T_wall_gas, T_gas, T_gas, P_tank, rho_gas, 0.021, 0.4, self.diam_in, "N2O") #relative to gas cv

        #NOTE: USE ambient properties for air, T_2 will be respective wall temperature (RK var)
        # (6) [natural convection] from atm to liq wall
        Q_dot_atm_to_liq_wall = solve_Q_dot_natural_convection_gas(self.T_atm, T_wall_liq, self.T_atm, self.P_atm, self.rho_atm, 0.59, 0.25, self.diam_out, "air") #relative to wall_liq cv
        # (7) [natural convection] from atm to gas wall
        Q_dot_atm_to_gas_wall = solve_Q_dot_natural_convection_gas(self.T_atm, T_wall_gas, self.T_atm, self.P_atm, self.rho_atm, 0.59, 0.25, self.diam_out, "air") #relative to wall_gas cv
        # (8) [conduction] from liq wall to gas wall 
        Q_dot_liq_wall_to_gas_wall = solve_Q_dot_conduction( (T_wall_liq-T_wall_gas), self.height_tank, self.k_w, self.diam_in, self.diam_out) #relative to wall_liq cv

        #(1) mass transfer from injector already solved
        m_dot_inj = spi_model(self.Cd_1, self.A_inj_1, P_tank, P_cc, rho_liq)
        #(2) mass transfer by condensation
        V_gas = m_gas/rho_gas
        m_dot_cond = solve_m_dot_condensed(T_gas, P_tank, V_gas, t)

        #(3) mass transfer by evaporation #NOTE: EMPIRICAL FACTOR E = 2.1E4 HERE TO CORRECTLY MODEL HEAT TRANSFER ON Q_DOT_LIQ_TO_SAT_SURF
        
        #NOTE: should be ok to just take the difference of the departure functions below
        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
        latent_heat_evap_l = preos_l.Hvap(T_liq)/MW
        h_liq = preos_l.H_dep_l/MW #+ n2o_ig.Cpg*(T_liq - T_REF)

        preos_sat = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_sat, P=P_tank)
        h_sat = preos_sat.H_dep_l/MW #+ n2o_ig.Cpg*(T_liq - T_REF)

        m_dot_evap = ((2.1E4)*Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) / (latent_heat_evap_l + (h_sat - h_liq) )
        #^ denom is constant, numerator fluxes like crazy!
        #print("m_dot_evap: ", ((2.1E4)*Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas), (2.1E4)*Q_dot_liq_to_sat_surf ,- Q_dot_sat_surf_to_gas)

        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
        latent_heat_cond_g = (-1)*preos_g.Hvap(T_gas)/MW

        m_dot_liq, m_dot_gas = solve_m_dot_liq_gas(m_dot_evap, m_dot_cond, m_dot_inj)
        #print(m_dot_evap, m_dot_cond, m_dot_inj)

        Q_dot_liq = Q_dot_liq_wall_to_liq -(2.1E4)*Q_dot_liq_to_sat_surf -m_dot_evap*latent_heat_evap_l
        Q_dot_gas = Q_dot_gas_wall_to_gas +Q_dot_sat_surf_to_gas -m_dot_cond*latent_heat_cond_g
    
        #print("Q_dots: ", Q_dot_liq, Q_dot_gas,  Q_dot_gas_wall_to_gas ,Q_dot_sat_surf_to_gas ,-m_dot_cond*latent_heat_cond_g)
        V_dot_liq = self.V_dot_liq_prev #initial guess for dV_dt_liq
        V_liq = self.V_tank - V_gas

        while np.abs(P_dot_error(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas)) > self.P_dot_err_tolerance:
            V_dot_liq = secant((lambda V_dot: P_dot_error(V_dot, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas )), V_dot_liq)
        
        ###solving wall nodes:
        height_dot = V_dot_liq / (0.25*np.pi*(self.diam_in**2))

        m_dot_liq_wall = self.rho_wall*(0.25*np.pi*height_dot*((self.diam_out**2)-(self.diam_in**2)))  #BUG: this might be a bit unstable w runge kutta steps?
        m_dot_gas_wall = -m_dot_liq_wall

        #print("m_dot_liq_wall should be - ", m_dot_liq_wall)


        #then solve the height of the gas wall
        h_gas_wall = V_gas / (0.25*np.pi*(self.diam_in**2))
        V_gas_wall = 0.25*np.pi*((self.diam_out**2)-(self.diam_in**2))*h_gas_wall
        m_gas_wall = self.rho_wall*V_gas_wall

        h_liq_wall = self.height_tank - h_gas_wall
        V_liq_wall = 0.25*np.pi*((self.diam_out**2)-(self.diam_in**2))*h_liq_wall
        m_liq_wall = self.rho_wall*V_liq_wall

        #NOTE: for T_dot_wall_liq:  IN: (6)      OUT: (4) AND (8)
        T_dot_wall_liq = ( Q_dot_atm_to_liq_wall - Q_dot_liq_wall_to_liq - Q_dot_liq_wall_to_gas_wall - m_dot_liq_wall*(0.15)*( T_wall_gas - T_wall_liq) ) / (0.15*m_liq_wall)

        #NOTE: for T_dot_wall_vap:  IN: (7) and (8)      OUT: (5)
        T_dot_wall_gas = ( Q_dot_atm_to_gas_wall - Q_dot_gas_wall_to_gas + Q_dot_liq_wall_to_gas_wall - m_dot_gas_wall*(0.15)*(T_wall_liq - T_wall_gas) ) / (0.15*m_gas_wall)

        #print("T_dot_wall_liq: ",  Q_dot_atm_to_liq_wall, - Q_dot_liq_wall_to_liq, - Q_dot_liq_wall_to_gas_wall, - m_dot_liq_wall*(0.15)*(T_wall_gas - T_wall_liq) )
        #print("T_dot_wall_gas: ",  Q_dot_atm_to_gas_wall, - Q_dot_gas_wall_to_gas, + Q_dot_liq_wall_to_gas_wall, - m_dot_gas_wall*(0.15)*(T_wall_liq - T_wall_gas) )

        #print("T_dot_wall_liq,gas: ", T_dot_wall_liq, T_dot_wall_gas)
        #print("wall masses! ", m_dot_gas_wall, m_dot_cond)

        print("\nlast single solve T_dot liq gas below:")
        T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas)


        print("sol derivs ", T_dot_liq, T_dot_gas, m_dot_liq, m_dot_gas, T_dot_wall_liq, T_dot_wall_gas)
        #print("heat transfer: ", Q_dot_liq, Q_dot_gas, (Q_dot_atm_to_liq_wall + Q_dot_liq_wall_to_liq - Q_dot_liq_wall_to_gas_wall + m_dot_liq_wall*(0.15)*( T_wall_gas - T_wall_liq)),(Q_dot_atm_to_gas_wall + Q_dot_gas_wall_to_gas + Q_dot_liq_wall_to_gas_wall + m_dot_gas_wall*(0.15)*(T_wall_liq - T_wall_gas)))
        #print("checking m_dot liq gas: ", m_dot_liq, m_dot_gas, m_dot_inj, m_dot_evap, m_dot_cond)

        return [T_dot_liq, T_dot_gas, m_dot_liq, m_dot_gas, T_dot_wall_liq, T_dot_wall_gas]

    def inst(self, P_cc):

        #TODO: check for liquid phase!!!!
        t = self.TIMESTEP
        y0 = [self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas]

        # NOTE: y is a vector, k1-k4 are derivatives of sol!
        k1 = self.system_of_liq_odes(t, y0, P_cc)

        y_2 = [y_i + self.TIMESTEP * k1_i / 2 for y_i, k1_i in zip(y0, k1)]
        k2 = self.system_of_liq_odes(t + self.TIMESTEP / 2, y_2, P_cc)

        y_3 = [y_i + self.TIMESTEP * k2_i / 2 for y_i, k2_i in zip(y0, k2)]
        k3 = self.system_of_liq_odes(t + self.TIMESTEP / 2, y_3, P_cc)

        y_4 = [y_i + self.TIMESTEP * k3_i for y_i, k3_i in zip(y0, k3)]
        k4 = self.system_of_liq_odes(t + self.TIMESTEP, y_4, P_cc)

        y =[ y_i + (self.TIMESTEP / 6) * (k1_i + 2*k2_i + 2*k3_i + k4_i) for y_i, k1_i,k2_i,k3_i,k4_i in zip(y0, k1,k2,k3,k4)]

        self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas = y

        #print("final y: ",self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas, "\n")
        print("final y: (expect +)", self.T_wall_liq - self.T_wall_gas, "\n")

        # (4) iteratively solve P_tank to update thermodynamic properties in each node
        #NOTE: this is just to update vals for downstream graphs
        self.rho_liq, self.rho_gas, self.P_tank = solve_thermo_params(self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.P_tank_prev, self.V_tank, self.volume_err_tolerance)

        #update stored vals for RK est and volumes
        self.P_tank_prev = self.P_tank
        self.V_dot_liq_prev = self.V_liq - self.m_liq/self.rho_liq
        self.V_liq = self.m_liq/self.rho_liq
        self.V_gas = self.V_tank - self.V_liq



t = 0
TIMESTEP = 1e-4

P_atm = 1e5 #Pa
T_atm = 286.5 #K
rho_atm = 1.225 #kg/m^3

m_nos = 20 #kg
P_tank = 45e5 #Pa
V_tank = 0.0354 #m^3

diam_out = 0.230 #m #NOTE: thesis didn't provide tank geometry, estimated based off of G type nos dimensions (approx equivalent nos mass to Karabeyoglu run tank)
diam_in = 0.215 #m
rho_wall = 2770 #kg/m^3
k_w = 237 #W/(m K)

Cd_1 = 0.425
A_inj_1 = 0.00003 #m^3 NOTE: GUESS
P_cc = 3e6 #1.03e6 + 68947.6 #Pa

inj_model = None #TODO: implement

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

P_sat_arr = []

T_liq_wall_arr = []
T_gas_wall_arr = []

#v_liq_arr = []
#v_vap_arr = []



###TODO: try solving single solve different ways!


try:
    while(t<1000*TIMESTEP):
        tank.inst(P_cc)
        t+=TIMESTEP 

        time_arr.append(t)
        #print("\n next timestep \n")
        
        P_tank_arr.append(tank.P_tank)
        m_tank_arr.append( (tank.m_liq+tank.m_gas) )
        m_liq_arr.append(tank.m_liq)
        m_gas_arr.append(tank.m_gas)



        T_liq_arr.append(tank.T_liq)
        T_gas_arr.append(tank.T_gas)

        T_liq_wall_arr.append(tank.T_wall_liq)
        T_gas_wall_arr.append(tank.T_wall_gas)

        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=tank.T_gas, P=tank.P_tank)
        T_sat = preos_g.Tsat(tank.P_tank)
        P_sat = preos_g.Psat(tank.T_gas)
        T_sat_arr.append(T_sat)
        P_sat_arr.append(P_sat)

except Exception as e:
    traceback.print_exc()


#print("Coolprop: ", CP.PhaseSI('T', 287.9284541887653, 'P', 4495139.921186801, "N2O"))

plt.subplot(1,4,1)
plt.scatter(time_arr,P_tank_arr,label = "tank")
plt.scatter(time_arr,P_sat_arr,label = "P_sat")
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Pa)')
plt.title('Pressure vs. Time')
plt.legend()
plt.grid(True)

plt.subplot(1,4,2)
plt.scatter(time_arr,m_tank_arr, label = "total mass")
plt.scatter(time_arr,m_liq_arr, label = "liquid")
plt.scatter(time_arr,m_gas_arr, label = "gas")
plt.xlabel('Time (s)')
plt.ylabel('Mass (kg)')
plt.title('Mass vs. Time')
plt.grid(True)

plt.subplot(1,4,3)
plt.scatter(time_arr,T_liq_arr, label = "liquid")
plt.scatter(time_arr,T_gas_arr, label = "gas")
plt.scatter(time_arr,T_sat_arr, label = "T_sat")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.title('Temperature vs. Time')
plt.legend()
plt.grid(True)

plt.subplot(1,4,4)
plt.scatter(time_arr,T_liq_wall_arr, label = "WALL liquid")
plt.scatter(time_arr,T_gas_wall_arr, label = "WALL gas")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.title('WALL Liquid Temperature vs. Time')
plt.legend()
plt.grid(True)

plt.show()


#variables at t=0 rk step 113:
#for ref V_dot_liq guess: -6.271942486202331e-09 
"""
m_liq = 18.405456101475963 
m_gas = 1.5932273922225426 
T_liq = 287.9284541887653
T_gas = 287.84453806778686
rho_liq = 784.2711376839605
rho_gas = 135.51112231223715
V_liq = 0.023468230887948742
V_gas = 0.01175717066641598 
P_tank = 4495139.921186801
m_dot_inj = -0.9399801595385376
m_dot_evap = 0.006204506778177737
m_dot_cond = 0
Q_dot_liq = 0.041625771932937 
Q_dot_gas = -0.028462026773206903
"""
"""
P_guess = 4500000.0
T_liq = 287.94493237008703
T_gas = 287.94493237008703
m_liq = 18.406207917835648
m_gas = 1.5937920821643532 
V_tank = 0.035225401554364724

P_guess = 4498592.6693985425 
T_liq = 287.9390159337461
T_gas =  287.91458387994976
m_liq = 18.405774049454863 
m_gas = 1.5937556870407577
V_tank =  0.035225401554364724

#print(P_guess, T_liq, T_gas, m_liq, m_gas, V_tank)
V_tank_error_arr = []
P_guess_arr = np.linspace(4500000.0*1.1, 4500000.0*.9, 100)
for x in P_guess_arr:
    y = V_tank_error(x, T_liq, T_gas, m_liq, m_gas, V_tank)
    #print(y)
    V_tank_error_arr = np.append(V_tank_error_arr , y)

plt.plot(P_guess_arr,V_tank_error_arr)
#plt.axhline(P_dot_err_tolerance, label = 'threshold of error to select solution to secant method', color = 'r')
plt.xlabel('P_guess')
plt.ylabel('V_tank_error')
plt.title('V_tank_error vs P_guess')
plt.grid(True)
#plt.legend()
plt.show()
"""










V_dot_arr = np.linspace(-1e-14, -1e-1, 100)
P_dot_error_arr = np.array([])

#NOTE: inputting V_dot_guess
"""
for x in V_dot_arr:
    y = P_dot_error(x, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, V_liq, V_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas)
    print(y)
    P_dot_error_arr = np.append(P_dot_error_arr, y)

plt.scatter(V_dot_arr,P_dot_error_arr)
#plt.axhline(P_dot_err_tolerance, label = 'threshold of error to select solution to secant method', color = 'r')
plt.xlabel('V_dot')
plt.ylabel('P_dot_error')
plt.title('P_dot_error vs V_dot')
plt.grid(True)
#plt.legend()
plt.show()
"""