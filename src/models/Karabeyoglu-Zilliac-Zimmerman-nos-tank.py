from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical

# Global Constants:
R_U = 8.31446 #J/(mol K) 

n2o = Chemical('nitrous oxide')

MW = (n2o.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o.omega - 0.26992*n2o.omega**2

TANK_DIAM = 0.0254*5.5 #m
CS_AREA = 0.25*np.pi*(TANK_DIAM**2) #m^2
g = 9.81 #m/s^2

#TODO: update with other injector model once we get this thing up
def spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit):
    m_dot_spi = (-1)*Cd_hem_spi_dyer * A_inj_ox * np.sqrt( 2 * rho_tank_exit * (P_1 - P_2)  )
    #print(rho_tank_exit)
    return m_dot_spi




#CoolProp missing some NOS models thermo has
def get_thermal_conductivity(T, P):
    n2o.T = T #K
    n2o.P = P #Pa
    return n2o.kl  #(W/(m K))

def get_viscosity(T, P): #(dynamic viscosity)
    n2o.T = T #K
    n2o.P = P #P
    return n2o.ViscosityLiquid.calculate_P(T,P, "LUCAS")  #(Pa s)





def solve_latent_heat_evap(T,P): 
    
    phase = CP.PropsSI('Phase', "T", T, "P", P, "N2O")
    phase_str = CP.PhaseSI("T", T, "P", P, "N2O")

    #print("phase: ", phase, phase_str, T, P)

    h_delta = (CP.PropsSI('H', 'T', T, 'Q', 1, "N2O") - CP.PropsSI('H', 'T', T, 'Q', 0, "N2O")) #J/kg

    if phase == 0: #subcooled fluid 
        h_evap = h_delta + (CP.PropsSI('H', 'T', T, 'Q', 0, "N2O") - CP.PropsSI('H', 'T', T, 'P', P, "N2O") ) #J/kg

    if phase == 6: #sat liq vapor
        h_evap = h_delta


    #print("PHASE: ", CP.PhaseSI("T", T, "P", P, "N2O"), "T", T, "P", P)

    print("latent heat: ", h_evap)
    return h_evap #J/kg


def solve_latent_heat_condens(T,P):
    return -solve_latent_heat_evap(T,P)



def solve_m_dot_condensed(T_gas, P, V_gas):

    P_sat = CP.PropsSI('P','T',T_gas,'Q',0,"N2O")

    if P > P_sat:
        return ((P-P_sat)*V_gas*MW)/( (R_U/MW)*T_gas*TIMESTEP)
        #NOTE: this might be derived from peng robinson eos and wrong

    else:
        return 0





def solve_Q_natural_convection(T_fluid, T_2, P, rho_1, c, n, fluid):

    k_1 = get_thermal_conductivity(T_fluid, P) #(W/(m K))
    Cp_1 = CP.PropsSI('C', 'T', T_fluid, 'P', P, fluid) #J/(kg K)
    visc_1 = get_viscosity(T_fluid, P) #(Pa s)

    dV_dT_P = n2o.VolumeLiquid.TP_dependent_property_derivative_T(T_fluid, P)  #TODO: CHECK THIS
    #print(f"\n dV_dT_P = {dV_dT_P}")
    beta = dV_dT_P*rho_1 

    Gr = ((TANK_DIAM**3)*(rho_1**2)*g*beta*np.abs(T_fluid - T_2) ) / (visc_1**2)
    Pr = (Cp_1*visc_1)/ k_1

    X = Gr*Pr
    h = c * (k_1/TANK_DIAM) * X**n

    Q_dot = h*CS_AREA*(T_fluid-T_2)

    return Q_dot


def solve_U_dot(T, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq , Q_dot_net, T_for_enthalpy_of_form):
    U_dot = m_dot_inj*CP.PropsSI('H','P',P_tank,'T',T,'N2O') - m_dot_evap*solve_latent_heat_evap(T_for_enthalpy_of_form,P_tank) + m_dot_cond*solve_latent_heat_condens(T_for_enthalpy_of_form,P_tank) - P_tank*V_dot_liq + Q_dot_net
    

    print("m_dot_evap!", m_dot_evap)
    print("inside solve U_dot", m_dot_inj*CP.PropsSI('H','P',P_tank,'T',T,'N2O') ,- m_dot_evap*solve_latent_heat_evap(T_for_enthalpy_of_form,P_tank) , + m_dot_cond*solve_latent_heat_condens(T_for_enthalpy_of_form,P_tank) ,- P_tank*V_dot_liq ,+ Q_dot_net)
    print("first term? ", m_dot_inj, CP.PropsSI('H','P',P_tank,'T',T,'N2O') )
    return U_dot


def solve_T_dot(t, T, rho, m, u_prev, U_dot, rho_dot, m_dot):

    partial_dP_drho_const_T = CP.PropsSI('d(U)/d(D)|T', 'T', T, 'D', rho, 'N2O') #TESTING
    Cv = CP.PropsSI('CVMASS', 'T', T, 'D', rho, 'N2O')

    #estimate u!
    u = u_prev + (t*U_dot)/m #NOTE: I THINK THIS MAKES SENSE BUT DOUBLE CHECK

    print("inside T_dot: ", (1/Cv)* (1/m) * (U_dot - u*m_dot) , U_dot , -u*m_dot)
    T_dot = (1/Cv)*( (1/m) * (U_dot - u*m_dot) - (partial_dP_drho_const_T* rho_dot) )
    return T_dot



def solve_Q_dot_conduction(delta_T, h_tank, k_w, diam_in, diam_out):

    L_w_cond = 0.5*h_tank
    Q_dot_conduction = k_w *(delta_T)*(0.25*np.pi*((diam_out**2)-(diam_in**2)))/L_w_cond
    
    return Q_dot_conduction


#NOTE: THIS WILL JUST RETURN PREV VALS IF t=0, which is good
def solve_thermo_properties(t, T_liq, T_gas, m_liq, m_gas, m_dot_liq, m_dot_gas, V_TANK, V_liq_prev, V_gas_prev, rho_liq_prev, rho_gas_prev, V_dot_liq_prev):

    
    # Initial guess for dV/dt_liq
    V_dot_liq_guess = V_dot_liq_prev  # BUG CHANGE THIS --> PASS IN PREVIOUS AS INPUT You can set this based on the previous timestep or system behavior 
    V_dot_liq_iter = -0.0000005
    
    volume_tolerance = 1e-4 
    pressure_tolerance = 2e2
    max_iterations = 4000 #TODO: define on max V_dot? can we thermodynamically figure out a maximum with some sort of hand calc assumption?

    for i in range(max_iterations):
        #print("V_DOT_LIQ_GUESS", V_dot_liq_guess)
        # Calculate dV/dt_vap based on dV/dt_liq
        V_dot_gas = -V_dot_liq_guess

        #use initial guess to solve future V_liq and V_gas with fwd euler method:
        #print("time: ", t, TIMESTEP)
        V_liq_est = V_liq_prev + t * V_dot_liq_guess
        V_gas_est = V_gas_prev + t * V_dot_gas
        #print("change in volume: ", V_liq_est, V_liq_prev, V_gas_est, V_gas_prev)
        
        # Calculate d_rho/dt for liquid and vapor
        rho_dot_liq = (1 / V_liq_est) * m_dot_liq - (m_liq / V_liq_est**2) * V_dot_liq_guess 
        rho_dot_gas = (1 / V_gas_est) * m_dot_gas - (m_gas / V_gas_est**2) * V_dot_gas    

        #print("checking assumption: ", (m_liq / V_liq_est**2),(m_gas / V_gas_est**2) )
        #print("masses",m_liq, m_gas)   
        
        # Update densities using Forward Euler method
        rho_liq = rho_liq_prev + t * rho_dot_liq
        rho_gas = rho_gas_prev + t * rho_dot_gas

        #print("checking dens change: ",rho_liq,rho_liq_prev, rho_gas,rho_gas_prev)

        #print("rho", (1 / V_gas_est) * m_dot_gas , - (m_gas / V_gas_est**2) * V_dot_gas ,rho_liq, rho_gas, t)
        
        # Calculate pressures using Peng-Robinson EOS #NOTE: try this first, if it doesnt look good try coolprop but that might not work not sure we can see 
        #pr_eos_liq = PR(T=T_liq, V=(MW/rho_liq), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
        P_liq = CP.PropsSI('P', 'T', T_liq, 'D', rho_liq, "N2O") #pr_eos_liq.P # 

        #pr_eos_gas = PR(T=T_gas, V=(MW/rho_gas), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
        P_gas = CP.PropsSI('P', 'T', T_gas, 'D', rho_gas, "N2O") #pr_eos_gas.P  

        #print("solved thermo properties 1 liq: ", rho_liq, V_liq_est, T_liq, P_liq, CP.PhaseSI("T",T_liq,"P",P_liq,"N2O"))
        #print("solved thermo properties 1 vap: ", rho_gas, V_gas_est, T_gas, P_gas, CP.PhaseSI("T",T_gas,"P",P_gas,"N2O"))
    
        
        # Check constraints:
        volume_constraint = abs(V_TANK - (m_liq / rho_liq + m_gas / rho_gas))
        pressure_constraint = abs(P_liq - P_gas)
        #print(m_liq , rho_liq , m_gas , rho_gas)
        #print(volume_constraint, pressure_constraint)
        #TODO: pressure derivative constraint?


        #print("vol const",volume_constraint, V_tank,- (m_liq / rho_liq + m_gas / rho_gas) )
        #print("pres const",pressure_constraint,P_liq , P_gas)
        
        # If both constraints are satisfied, return the solution
        if volume_constraint < volume_tolerance and pressure_constraint < pressure_tolerance:
            #print("number of loops: ", i)
            return P_gas, rho_liq, rho_gas, V_liq_est, V_gas_est, V_dot_liq_guess, V_dot_gas, rho_dot_liq, rho_dot_gas #, V_dot_liq_guess

        V_dot_liq_guess += V_dot_liq_iter
        i+=1

    #print("solved thermo properties (fail) liq: ", rho_liq, V_liq_est, T_liq, P_tank, CP.PhaseSI("T",T_liq,"P",P_liq,"N2O"))
    #print("solved thermo properties (fail) vap: ", rho_gas, V_gas_est, T_gas, P_tank, CP.PhaseSI("T",T_gas,"P",P_gas,"N2O"))

    print("fail to conv", P_gas, P_liq, rho_liq, rho_gas, V_liq_est, V_gas_est, V_dot_liq_guess, V_dot_gas, rho_dot_liq, rho_dot_gas)
    #print(volume_constraint, pressure_constraint)
    raise Exception("Failed to converge")


class LiquidTankODES:
    def __init__(self, m_dot_liq, m_dot_gas):
        # Initialize previous values
        self.m_dot_liq_prev = m_dot_liq
        self.m_dot_gas_prev = m_dot_gas

        self.V_dot_liq_prev = 0

        self.thermo_sol = []


    def system_of_liq_odes(self, t, y, constants):

        T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas = y  # Unpack state variables
        print("unpacked state var: ", T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas)

        V_TANK, V_liq_prev, V_gas_prev, rho_liq_prev, rho_gas_prev, rho_wall, diam_out, diam_in, u_prev_liq, u_prev_gas, k_w, height_tank, T_atm, P_atm, rho_atm = constants

        #print(m_liq, m_gas)
        #solve thermo properties?

        P_tank, rho_liq, rho_gas, V_liq, V_gas, V_dot_liq, V_dot_gas, rho_dot_liq, rho_dot_gas = solve_thermo_properties(t, T_liq, T_gas, m_liq, m_gas, self.m_dot_liq_prev, self.m_dot_gas_prev, V_TANK, V_liq_prev, V_gas_prev, rho_liq_prev, rho_gas_prev, self.V_dot_liq_prev)
    
        print("solved thermo prop: ", P_tank, rho_liq, rho_gas,)
        self.V_dot_liq_prev = V_dot_liq
        #print("solved thermo properties liq: ", rho_liq, V_liq, T_liq, P_tank, CP.PhaseSI("T",T_liq,"P",P_tank,"N2O"))
        #print("solved thermo properties vap: ", rho_gas, V_gas, T_gas, P_tank, CP.PhaseSI("T",T_gas,"P",P_tank,"N2O"))
        
        self.thermo_sol= [P_tank, rho_liq, rho_gas, V_liq, V_gas, V_dot_liq, V_dot_gas, rho_dot_liq, rho_dot_gas]

        #NOTE: BUG BUG heat transfer signs likely wrong, I think everything should be + because Q_dot eqns will return signed vals

        #solve heat transfer terms
        T_sat = CP.PropsSI('T','P',P_tank,'Q',0,'N2O')

        # (2) from saturated surface to gas
        Q_dot_sat_surf_to_gas = solve_Q_natural_convection(T_gas, T_sat, P_tank, rho_gas, 0.15, 0.333, "N2O") #relative to gas cv
        # (3)  from liq to saturated surface
        Q_dot_liq_to_sat_surf = solve_Q_natural_convection(T_liq, T_sat, P_tank, rho_liq, 0.15, 0.333, "N2O") #relative to liq cv
        # (4) [natural convection] from liq wall to liq
        Q_dot_liq_wall_to_liq = solve_Q_natural_convection(T_liq, T_wall_liq, P_tank, rho_liq, 0.59, 0.25, "N2O") #relative to liq wall cv
        # (5) [natural convection] from gas wall to gas
        Q_dot_gas_wall_to_gas = solve_Q_natural_convection(T_gas, T_wall_gas, P_tank, rho_gas, 0.59, 0.25, "N2O") #relative to gas wall cv

        #print( Q_dot_sat_surf_to_gas, Q_dot_liq_to_sat_surf, Q_dot_liq_wall_to_liq, Q_dot_gas_wall_to_gas)

        #NOTE: USE ambient properties for air, T_2 will be respective wall temperature (RK var)
        # (6) [natural convection] from atm to liq wall
        Q_dot_atm_to_liq_wall = solve_Q_natural_convection(T_atm, T_wall_liq, P_atm, rho_atm, 0.59, 0.25, "air") #relative
        # (7) [natural convection] from atm to vap wall
        Q_dot_atm_to_gas_wall = solve_Q_natural_convection(T_atm, T_wall_gas, P_atm, rho_atm, 0.59, 0.25, "air")
        # (8) [conduction] from liq wall to vap wall
        Q_dot_liq_wall_to_gas_wall = solve_Q_dot_conduction( (T_wall_liq-T_wall_gas), height_tank, k_w, diam_in, diam_out) #TODO: check sign convention

        #NOTE: all above start at 0 since starting at thermal eq assumption
        #print( Q_dot_atm_to_liq_wall, Q_dot_atm_to_gas_wall, Q_dot_liq_wall_to_gas_wall)

        #solve mass transfer NOTE: all three values sign convention relative to liquid cv
        #(1) mass transfer from injector already solved
        m_dot_inj = spi_model(Cd_1, A_inj_1, P_tank, P_cc, rho_liq)
        #(2) mass transfer by condensation
        m_dot_cond = solve_m_dot_condensed(T_gas, P_tank, V_gas)
        #(3) mass transfer by evaporation
        print("checking temps? ", T_liq, T_gas)
        m_dot_evap = (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) / (solve_latent_heat_evap(T_liq, P_tank)) #NOTE: P sure its T_liq but might be wrong


        m_dot_liq = m_dot_evap + m_dot_cond + m_dot_inj
        m_dot_gas = (-1)* (m_dot_evap - m_dot_cond) #convert sign convention from liq cv to gas cv

        #print("net mass flow rates l,g",m_dot_liq, m_dot_gas)

        self.m_dot_liq_prev = m_dot_liq
        self.m_dot_gas_prev = m_dot_gas

        #TODO: solve m_dot_wall
        delta_height = (4*(V_liq - V_liq_prev)/(np.pi*(diam_in**2)))
        m_dot_wall = 0
        if t != 0:
            m_dot_wall = rho_wall*(0.25*np.pi*delta_height*((diam_out**2)-(diam_in**2)))/t

        #BUG: this is likely wrong ^

        #update m_wall liq and vap
        m_wall_liq = 0.25*np.pi*((diam_out**2)-(diam_in**2)) * ( (4*V_liq)/(np.pi*diam_in**2) ) * rho_wall
        m_wall_gas = 0.25*np.pi*((diam_out**2)-(diam_in**2)) * ( (4*V_gas)/(np.pi*diam_in**2) ) * rho_wall
        

        #NOTE: for T_dot_wall_liq:  IN: (6)      OUT: (4) AND (8)
        T_dot_wall_liq = ( Q_dot_atm_to_liq_wall - Q_dot_liq_wall_to_liq - Q_dot_liq_wall_to_gas_wall +  m_dot_wall*(0.15)*( T_wall_gas - T_wall_liq) ) / (0.15*m_wall_liq)

        #NOTE: for T_dot_wall_vap:  IN: (7) and (8)      OUT: (5)
        T_dot_wall_gas = ( Q_dot_atm_to_gas_wall - Q_dot_gas_wall_to_gas - Q_dot_liq_wall_to_gas_wall +  m_dot_wall*(0.15)*(T_wall_liq - T_wall_gas) ) / (0.15*m_wall_gas)
        
        #BUG BUG BUG ^^^ CHECK THAT ONCE PAST FIRST TIMESTEP

        #NOTE: for Q_dot_liq:   IN:(2)       OUT: (3)
        U_dot_liq = solve_U_dot(T_liq, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq , (Q_dot_liq_wall_to_liq - Q_dot_liq_to_sat_surf), T_liq )

        #NOTE: for Q_dot_liq:   IN:(5)       OUT: (NONE)
        U_dot_gas = solve_U_dot(T_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_gas, Q_dot_gas_wall_to_gas, T_liq)


        T_dot_liq = solve_T_dot(t, T_liq, rho_liq, m_liq, u_prev_liq, U_dot_liq, rho_dot_liq, m_dot_liq)


        T_dot_gas = solve_T_dot(t, T_gas, rho_gas, m_gas, u_prev_gas, U_dot_gas, rho_dot_gas, m_dot_gas)

        print("Q dot: ", Q_dot_gas_wall_to_gas, (Q_dot_liq_wall_to_liq - Q_dot_liq_to_sat_surf))
        print("U dot: ",U_dot_liq, U_dot_gas )
        print("out derivs: ",T_dot_liq, T_dot_gas, m_dot_liq, m_dot_gas, T_dot_wall_liq, T_dot_wall_gas)

        return [T_dot_liq, T_dot_gas, m_dot_liq, m_dot_gas, T_dot_wall_liq, T_dot_wall_gas]








class model():

    #self.V_tank, self.V_liq, self.V_vap, self.rho_liq, self.rho_vap, self.rho_wall, self.diam_out, self.diam_in, self.u_liq, self.u_vap, self.k_w, ((4*self.V_tank)/(np.pi*self.diam_in**2))
    def __init__(self,TIMESTEP, m_nos, P_tank, V_tank, diam_out, diam_in, rho_wall, k_w, Cd_1, A_inj_1, P_cc, T_atm, P_atm, rho_atm, inj_model):
        
        #tank geometry
        self.rho_wall = rho_wall
        self.k_w = k_w
        self.diam_out = diam_out
        self.diam_in = diam_in

        #assuming tank wall starts in thermal equilibrium with environment?
        self.T_wall_liq = T_atm
        self.T_wall_vap = T_atm

        self.T_atm = T_atm
        self.P_atm = P_atm
        self.rho_atm = rho_atm

        #injector input
        self.Cd_1 = Cd_1
        self.P_cc = P_cc
        self.A_inj_1 = A_inj_1

        self.P_tank = P_tank

        #NOTE: Setup - Assuming Tank Starts at Sat Conditions
        self.TIMESTEP = TIMESTEP

        self.V_tank = V_tank
        self.m_nos = m_nos

        rho_sat_vap = CP.PropsSI('D', 'P', self.P_tank, 'Q', 1, "N2O")
        rho_sat_liq = CP.PropsSI('D', 'P', self.P_tank, 'Q', 0, "N2O")
        rho_sat_tank = self.m_nos/V_tank

        x_tank = ( (1/rho_sat_tank)-(1/rho_sat_liq) ) / ( (1/rho_sat_vap)-(1/rho_sat_liq) )

        ### initial cv values from initial sat liq vap state

        #vapor cv
        self.m_vap = x_tank*self.m_nos

        self.T_vap = CP.PropsSI('T', 'P', self.P_tank, 'Q', 1, "N2O") + 0.005 #Perturb to start close to equillibrium
        #print("PHASE CHECK, ",CP.PhaseSI("T", self.T_vap, "P", self.P_tank, "N2O"))

        self.rho_vap = CP.PropsSI('D','P',self.P_tank,'T',self.T_vap,"N2O")
        self.v_vap = 1/self.rho_vap
        self.V_vap = self.v_vap*self.m_vap
        self.u_vap = CP.PropsSI('U','P',self.P_tank,'T',self.T_vap,"N2O")


        #liquid cv
        self.m_liq = self.m_nos-self.m_vap

        self.T_liq = CP.PropsSI('T', 'P', self.P_tank, 'Q', 0, "N2O") - 0.005 #Perturb to start close to equillibrium

        self.rho_liq = CP.PropsSI('D','P',self.P_tank,'T',self.T_liq,"N2O")
        self.v_liq = 1/self.rho_liq
        self.V_liq = self.V_tank-self.V_vap
        self.u_liq = CP.PropsSI('U','P',self.P_tank,'T',self.T_liq,"N2O")

                
        self.P_vap_prev = P_tank
        self.P_liq_prev = P_tank

        self.rho_vap_prev = (1/self.v_vap)
        self.rho_liq_prev = (1/self.v_liq)

        self.rho_exit = 1/self.v_liq

        self.t = 0

        self.m_dot_liq = 0
        self.m_dot_vap = 0

        #print(self.V_liq, self.T_liq, self.rho_liq, self.P_tank)
        

        

        



    def inst(self, P_cc):

        self.t+=self.TIMESTEP
        
        if self.m_liq >= 0: #liquid phase!
            
            #enter 6 ode system

            

            constants = [ self.V_tank, self.V_liq, self.V_vap, self.rho_liq, self.rho_vap, self.rho_wall, self.diam_out, self.diam_in, self.u_liq, self.u_vap, self.k_w, ((4*self.V_tank)/(np.pi*self.diam_in**2)), self.T_atm, self.P_atm, self.rho_atm ]
            y0 =  [self.T_liq, self.T_vap, self.m_liq, self.m_vap, self.T_wall_liq, self.T_wall_vap ]

            #print(self.P_tank, self.V_tank, self.V_liq, self.V_vap, self.rho_liq, self.rho_vap, self.rho_wall, self.diam_out, self.diam_in, self.u_liq, self.u_vap, self.k_w, ((4*self.V_tank)/(np.pi*self.diam_in**2)), self.T_atm, self.P_atm, self.rho_atm )

            #setup odes:
            odes = LiquidTankODES(self.m_dot_liq, self.m_dot_vap)
            solution = solve_ivp(odes.system_of_liq_odes, [0, self.TIMESTEP], y0, args=(constants,), method='RK45', rtol=1e-12, atol=1e-12) #NOTE: THIS IS PROBABLY TOO SLOW, I SET THIS HIGH TOLERANCE FOR THE LAST TRY
            

            # Extract the final state from the solution
            final_state = solution.y[:, -1]  # Get the last column, which represents the state at the final time step

            # Unpack the final state into your variables
            self.T_liq, self.T_gas, self.m_liq, self.m_vap, self.T_wall_liq, self.T_wall_gas = final_state


            #TODO: solve thermo properties
            self.P_tank, self.rho_liq, self.rho_gas, self.V_liq, self.V_gas, V_dot_liq, V_dot_gas, rho_dot_liq, rho_dot_gas = odes.thermo_sol
            
            #print("checking volume and sum", self.V_liq, self.V_gas, self.V_liq+ self.V_gas)
            #solve thermo properties or return them from ^ (latter prefered)
            #update tank properties!


            
            


        else: #vapor phase

            return 0










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
tank = model(TIMESTEP, m_nos, P_tank, V_tank, diam_out, diam_in, rho_wall, k_w, Cd_1, A_inj_1, P_cc, T_atm, P_atm, rho_atm, inj_model)
while(t<5*TIMESTEP):
    tank.inst(P_cc)
    t+=TIMESTEP
    print(t, tank.P_tank, tank.V_liq,"\n")
    #print("\n")