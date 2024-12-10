from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical
import matplotlib.pyplot as plt

# Global Constants:
R_U = 8.31446 #J/(mol K) 

n2o = Chemical('nitrous oxide')

MW = (n2o.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o.omega - 0.26992*n2o.omega**2

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
        kk = kk + 1
    x = x2
    return x

#CoolProp missing some NOS models thermo has
def get_thermal_conductivity(T, P):
    n2o.T = T #K
    n2o.P = P #Pa
    return n2o.kl  #(W/(m K))

def get_viscosity(T, P): #(dynamic viscosity)
    n2o.T = T #K
    n2o.P = P #P
    return n2o.ViscosityLiquid.calculate_P(T,P, "LUCAS")  #(Pa s)

def solve_Q_dot_natural_convection(T_1, T_2, P, rho_2, c, n, tank_diam, g, fluid):

    k_1 = get_thermal_conductivity(T_2, P) #(W/(m K))
    Cp_1 = CP.PropsSI('C', 'T', T_2, 'P', P, fluid) #J/(kg K)
    visc_1 = get_viscosity(T_2, P) #(Pa s)

    dV_dT_P = n2o.VolumeLiquid.TP_dependent_property_derivative_T(T_2, P)  #TODO: CHECK THIS
    #print(f"\n dV_dT_P = {dV_dT_P}")
    beta = dV_dT_P*rho_2 

    Gr = ((tank_diam**3)*(rho_2**2)*g*beta*np.abs(T_1 - T_2) ) / (visc_1**2)
    Pr = (Cp_1*visc_1)/ k_1

    X = Gr*Pr
    h = c * (k_1/tank_diam) * X**n

    Q_dot = h*(0.25*np.pi*(tank_diam**2))*(T_1-T_2)

    return Q_dot

def solve_Q_dot_conduction(delta_T, h_tank, k_w, diam_in, diam_out):

    L_w_cond = 0.5*h_tank
    Q_dot_conduction = k_w *(delta_T)*(0.25*np.pi*((diam_out**2)-(diam_in**2)))/L_w_cond
    
    return Q_dot_conduction

#TODO: update with other injector model once we get this thing up
def spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit):
    m_dot_spi = (-1)*Cd_hem_spi_dyer * A_inj_ox * np.sqrt( 2 * rho_tank_exit * (P_1 - P_2)  )
    #print(rho_tank_exit)
    return m_dot_spi

def solve_m_dot_condensed(T_gas, P, V_gas, t):

    P_sat = CP.PropsSI('P','T',T_gas,'Q',0,"N2O")

    if P > P_sat:
        return ((P-P_sat)*V_gas*MW)/( (R_U/MW)*T_gas*t)
        #NOTE: this might be derived from peng robinson eos and wrong

    else:
        return 0

def solve_latent_heat_evap(T,P): 
    
    phase = CP.PropsSI('Phase', "T", T, "P", P, "N2O")

    h_delta = (CP.PropsSI('H', 'T', T, 'Q', 1, "N2O") - CP.PropsSI('H', 'T', T, 'Q', 0, "N2O")) #J/kg

    if phase == 0: #subcooled fluid 
        h_evap = h_delta + (CP.PropsSI('H', 'T', T, 'Q', 0, "N2O") - CP.PropsSI('H', 'T', T, 'P', P, "N2O") ) #J/kg

    if phase == 6: #sat liq gasor
        h_evap = h_delta

    return h_evap #J/kg

def solve_U_dot(T, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq , Q_dot_net, T_for_enthalpy_of_form):
    U_dot = m_dot_inj*CP.PropsSI('H','P',P_tank,'T',T,'N2O') + m_dot_evap*solve_latent_heat_evap(T_for_enthalpy_of_form,P_tank) + m_dot_cond*(-1)*solve_latent_heat_evap(T_for_enthalpy_of_form,P_tank) - P_tank*V_dot_liq + Q_dot_net #NOTE: (-1)* to convert latent heat evap to condensation

    return U_dot

def V_tank_error(P_guess, T_liq, T_gas, m_liq, m_gas, V_tank):

    rho_liq = CP.PropsSI('D', 'T', T_liq, 'P', P_guess, 'N2O')
    rho_gas = CP.PropsSI('D', 'T', T_gas, 'P', P_guess, 'N2O')

    V_tank_est = (m_liq/rho_liq) + (m_gas/rho_gas)

    return V_tank_est - V_tank 

def solve_thermo_params(T_liq, T_gas, m_liq, m_gas, P_tank_prev, volume_err_tolerance):

    P_tank = P_tank_prev #initial guess for pressure
    while np.abs(V_tank_error(P_tank, T_liq, T_gas, m_liq, m_gas) ) > volume_err_tolerance:
        P_tank = secant((lambda P: V_tank_error(P, T_liq, T_gas, m_liq, m_gas)), P_tank)

    rho_liq = CP.PropsSI('D', 'T', T_liq, 'P', P_tank, 'N2O')
    rho_gas = CP.PropsSI('D', 'T', T_gas, 'P', P_tank, 'N2O')

    return rho_liq, rho_gas, P_tank

def single_solve_T_dot_liq_gas(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas):

    m_dot_liq = m_dot_evap + m_dot_cond + m_dot_inj
    m_dot_gas = (-1)*(m_dot_evap - m_dot_cond) #convert sign convention from liq cv to gas cv

    V_dot_gas = -V_dot_liq

    V_liq = m_liq/rho_liq
    V_gas = m_gas/rho_gas

    rho_dot_liq = (1 / V_liq) * m_dot_liq - (m_liq / V_liq**2) *  V_dot_liq
    rho_dot_gas = (1 / V_gas) * m_dot_gas - (m_gas / V_gas**2) * V_dot_gas  

    U_dot_liq = solve_U_dot(T_liq, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq , Q_dot_liq, T_liq )
    U_dot_gas = solve_U_dot(T_gas, P_tank, 0, (-1)*m_dot_evap, m_dot_cond, V_dot_gas, Q_dot_gas, T_liq)


    partial_du_drho_const_T_liq = CP.PropsSI('d(U)/d(D)|T', 'T', T_liq, 'D', rho_liq, 'N2O') #NOTE: TESTING, not sure if cp can do this
    Cv_liq = CP.PropsSI('CVMASS', 'T', T_liq, 'D', rho_liq, 'N2O')
    u_liq = CP.PropsSI('U', 'T', T_liq, 'D', rho_liq, 'N2O')

    partial_du_drho_const_T_gas = CP.PropsSI('d(U)/d(D)|T', 'T', T_gas, 'D', rho_gas, 'N2O') #NOTE: TESTING, not sure if cp can do this
    Cv_gas = CP.PropsSI('CVMASS', 'T', T_gas, 'D', rho_gas, 'N2O')
    u_gas = CP.PropsSI('U', 'T', T_gas, 'D', rho_gas, 'N2O')
    
    T_dot_liq = (1/Cv_liq)*( (1/m_liq) * (U_dot_liq - u_liq*m_dot_liq) - (partial_du_drho_const_T_liq* rho_dot_liq) )
    T_dot_gas = (1/Cv_gas)*( (1/m_gas) * (U_dot_gas - u_gas*m_dot_gas) - (partial_du_drho_const_T_gas* rho_dot_gas) )
    return T_dot_liq, T_dot_gas

def P_dot_error(V_dot_guess, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas):    

    m_dot_liq = m_dot_evap + m_dot_cond + m_dot_inj
    m_dot_gas = (-1)*(m_dot_evap - m_dot_cond) #convert sign convention from liq cv to gas cv

    V_liq = m_liq/rho_liq
    V_gas = m_gas/rho_gas

    T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_guess, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas)

    rho_dot_liq = (1 / V_liq) * m_dot_liq - (m_liq / V_liq**2) *  V_dot_guess
    rho_dot_gas = (1 / V_gas) * m_dot_gas - (m_gas / V_gas**2) * V_dot_guess (-1)  

    partial_dP_dT_const_rho_liq = CP.PropsSI('d(P)/d(T)|D', 'T', T_liq, 'D', rho_liq, 'N2O')
    partial_dP_drho_const_T_liq = CP.PropsSI('d(P)/d(D)|T', 'T', T_liq, 'D', rho_liq, 'N2O')

    P_dot_liq = partial_dP_dT_const_rho_liq*T_dot_liq + partial_dP_drho_const_T_liq*rho_dot_liq

    partial_dP_dT_const_rho_gas = CP.PropsSI('d(P)/d(T)|D', 'T', T_gas, 'D', rho_gas, 'N2O')
    partial_dP_drho_const_T_gas = CP.PropsSI('d(P)/d(D)|T', 'T', T_gas, 'D', rho_gas, 'N2O')

    P_dot_gas = partial_dP_dT_const_rho_gas*T_dot_gas + partial_dP_drho_const_T_gas*rho_dot_gas


    return P_dot_liq - P_dot_gas

def iteratively_solve_T_dot_liq_and_gas(m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas, V_dot_liq_prev, P_dot_err_tolerance):

    V_dot_liq = V_dot_liq_prev #initial guess for dV_dt_liq

    while np.abs(P_dot_error(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas) ) > P_dot_err_tolerance:
        V_dot_liq = secant((lambda V_dot: P_dot_error(V_dot, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas )), V_dot_liq)

    T_dot_liq, T_dot_gas = single_solve_T_dot_liq_gas(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas)

    return T_dot_liq, T_dot_gas


def system_of_liq_odes(t, y, constants):

    T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas = y  # Unpack state variables

    T_atm, P_atm, rho_atm, height_tank, k_w, diam_in, diam_out, Cd_1, A_inj_1, P_cc, V_tank, rho_wall, P_tank_prev, V_dot_liq_prev, volume_err_tolerance, P_dot_err_tolerance = constants

    #solve thermodynamic properties with method of (4), since T, m will change so will the dependent thermodynamic properties!!!
    rho_liq, rho_gas, P_tank = solve_thermo_params(T_liq, T_gas, m_liq, m_gas, P_tank_prev, volume_err_tolerance)

    #solve heat transfer terms
    T_sat = CP.PropsSI('T','P',P_tank,'Q',0,'N2O')

    #NOTE: All heat transfer terms solved by solve_Q_dot_natural_convection have signs relative to the fluid cv not the surface

    # (2) from saturated surface to gas
    #TODO: FIX NAMING ITS Q DOT 
    Q_dot_sat_surf_to_gas = solve_Q_dot_natural_convection(T_sat, T_gas, P_tank, rho_gas, 0.15, 0.333, "N2O") #relative to gas cv
    # (3)  from liq to saturated surface
    Q_dot_liq_to_sat_surf = solve_Q_dot_natural_convection(T_sat, T_liq, P_tank, rho_liq, 0.15, 0.333, "N2O") #relative to liq cv
    # (4) [natural convection] from liq wall to liq
    Q_dot_liq_wall_to_liq = solve_Q_dot_natural_convection(T_wall_liq, T_liq, P_tank, rho_liq, 0.59, 0.25, "N2O") #relative to liq wall cv       
    # (5) [natural convection] from gas wall to gas
    Q_dot_gas_wall_to_gas = solve_Q_dot_natural_convection(T_wall_gas, T_gas, P_tank, rho_gas, 0.59, 0.25, "N2O") #relative to gas wall cv

    #print( Q_dot_sat_surf_to_gas, Q_dot_liq_to_sat_surf, Q_dot_liq_wall_to_liq, Q_dot_gas_wall_to_gas)

    #NOTE: USE ambient properties for air, T_2 will be respective wall temperature (RK var)
    # (6) [natural convection] from atm to liq wall
    Q_dot_atm_to_liq_wall = solve_Q_dot_natural_convection(T_wall_liq, T_atm, P_atm, rho_atm, 0.59, 0.25, "air") #relative
    # (7) [natural convection] from atm to gas wall
    Q_dot_atm_to_gas_wall = solve_Q_dot_natural_convection(T_wall_gas, T_atm, P_atm, rho_atm, 0.59, 0.25, "air")
    # (8) [conduction] from liq wall to gas wall
    Q_dot_liq_wall_to_gas_wall = solve_Q_dot_conduction( (T_wall_gas-T_wall_liq), height_tank, k_w, diam_in, diam_out) #relative to liquid cv


    #(1) mass transfer from injector already solved
    m_dot_inj = spi_model(Cd_1, A_inj_1, P_tank, P_cc, rho_liq)
    #(2) mass transfer by condensation
    V_gas = m_gas*rho_gas
    m_dot_cond = solve_m_dot_condensed(T_gas, P_tank, V_gas, t)
    #(3) mass transfer by evaporation
    m_dot_evap = (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) / (solve_latent_heat_evap(T_liq, P_tank)) #NOTE: P sure its T_liq but might be wrong

    m_dot_liq = m_dot_evap + m_dot_cond + m_dot_inj
    m_dot_gas = (-1)*(m_dot_evap - m_dot_cond) #convert sign convention from liq cv to gas cv

    ####solve wall properties!!! then heat transfer terms

    #TODO: solve m_dot_wall

    Q_dot_liq = Q_dot_liq_wall_to_liq - Q_dot_liq_to_sat_surf
    Q_dot_gas = Q_dot_gas_wall_to_gas

    V_dot_liq = V_dot_liq_prev #initial guess for dV_dt_liq

    while np.abs(P_dot_error(V_dot_liq, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas) ) > P_dot_err_tolerance:
        V_dot_liq = secant((lambda V_dot: P_dot_error(V_dot, m_liq, m_gas, T_liq, T_gas, rho_liq, rho_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, Q_dot_liq, Q_dot_gas )), V_dot_liq)

    delta_height = V_dot_liq / (0.25*np.pi*(diam_in**2))
    m_dot_liq_wall = rho_wall*(0.25*np.pi*delta_height*((diam_out**2)-(diam_in**2)))/t
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

    def __init__(self, m_nos, P_tank, P_cc, P_atm, T_atm, rho_atm, V_tank, diam_out, diam_in, rho_wall, k_w, Cd_1, A_inj_1, volume_err_tolerance, P_dot_err_tolerance):
        
        self.volume_err_tolerance = volume_err_tolerance
        self.P_dot_err_tolerance = P_dot_err_tolerance

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

        self.T_wall_liq = T_atm
        self.T_wall_gas = T_atm

        #injector input
        self.Cd_1 = Cd_1
        self.P_cc = P_cc
        self.A_inj_1 = A_inj_1

        self.P_tank = P_tank

        #NOTE: Setup Thermo Properties - Assuming Tank Starts NEAR Sat Conditions

        self.V_tank = V_tank
        self.m_nos = m_nos

        rho_sat_gas = CP.PropsSI('D', 'P', self.P_tank, 'Q', 1, "N2O")
        rho_sat_liq = CP.PropsSI('D', 'P', self.P_tank, 'Q', 0, "N2O")
        rho_sat_tank = self.m_nos/V_tank

        x_tank = ( (1/rho_sat_tank)-(1/rho_sat_liq) ) / ( (1/rho_sat_gas)-(1/rho_sat_liq) )

        #gas cv
        self.m_gas = x_tank*self.m_nos

        self.T_gas = CP.PropsSI('T', 'P', self.P_tank, 'Q', 1, "N2O") + 0.1 #Perturb to start close to equillibrium but as a gasor

        self.rho_gas = CP.PropsSI('D','P',self.P_tank,'T',self.T_gas,"N2O")
        v_gas = 1/self.rho_gas
        self.V_gas = v_gas*self.m_gas

        #liquid cv
        self.m_liq = self.m_nos-self.m_gas

        self.T_liq = CP.PropsSI('T', 'P', self.P_tank, 'Q', 0, "N2O") - 0.1 #Perturb to start close to equillibrium

        self.rho_liq = CP.PropsSI('D','P',self.P_tank,'T',self.T_liq,"N2O")
        v_liq = 1/self.rho_liq
        self.V_liq = v_liq*self.m_liq

        #self.rho_exit = 1/v_liq

        self.V_tank = self.V_liq+self.V_gas # "what are you going to do if the aluminum is too small? water it? give it sunlight? let it grow?"
        self.height_tank = self.V_tank/(0.25*np.pi*(diam_in**2))


    def inst(self, t, y, constants):
        
        #TODO: check for liquid phase!!!!

        # (1) enter Explicit Runge Kutta 6 ode system
        #T_atm, P_atm, rho_atm, height_tank, k_w, diam_in, diam_out, Cd_1, A_inj_1, P_cc, V_tank, rho_wall, P_tank_prev, V_dot_liq_prev, volume_err_tolerance, P_dot_err_tolerance
        constants = [self.T_atm, self.P_atm, self.rho_atm, self.height_tank, self.k_w, self.diam_in, self.diam_out, self.Cd_1, self.A_inj_1, self.P_cc, self.V_tank, self.rho_wall, self.P_tank_prev, self.V_dot_liq_prev, self.volume_err_tolerance, self.P_dot_err_tolerance]
        y0 =  [self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas ]

        #solve with scipy
        solution = solve_ivp(system_of_liq_odes, [0, self.TIMESTEP], y0, args=(constants,), method='RK45', rtol=1e-12, atol=1e-12) #NOTE: THIS IS PROBABLY TOO SLOW, I SET THIS HIGH TOLERANCE FOR THE LAST TRY
            

        # Extract the final state from the solution
        final_state = solution.y[:, -1]  # Get the last column, which represents the state at the final time step

        # Unpack the final state into your variables
        self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.T_wall_liq, self.T_wall_gas = final_state

        # (4) iteratively solve P_tank to update thermodynamic properties in each node
        self.rho_liq, self.rho_gas, self.P_tank = solve_thermo_params(self.T_liq, self.T_gas, self.m_liq, self.m_gas, self.P_tank_prev, self.volume_err_tolerance)

        #update stored vals for RK est and volumes
        self.P_tank_prev = self.P_tank
        self.V_dot_liq_prev = self.V_liq - self.m_liq/self.rho_liq
        self.V_liq = self.m_liq/self.rho_liq
        self.V_gas = self.V_tank - self.V_liq
