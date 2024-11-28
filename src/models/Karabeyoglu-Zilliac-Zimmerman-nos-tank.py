from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical

# Global Constants:
P_ATM = 1e5 #Pa
T_ATM = 273.15 + 15 #K

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

    #print("phase: ", phase)

    h_delta = (CP.PropsSI('H', 'T', T, 'Q', 1, "N2O") - CP.PropsSI('H', 'T', T, 'Q', 0, "N2O")) #J/kg

    if phase == 0: #subcooled fluid 
        h_evap = h_delta + (CP.PropsSI('H', 'T', T, 'Q', 0, "N2O") - CP.PropsSI('H', 'T', T, 'P', P, "N2O") ) #J/kg

    if phase == 6: #sat liq vapor
        h_evap = h_delta


    #print("PHASE: ", CP.PhaseSI("T", T, "P", P, "N2O"), "T", T, "P", P)


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


def solve_U_dot(T, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq , Q_dot_net):
    U_dot = m_dot_inj*CP.PropsSI('H','P',P_tank,'T',T,'N2O') - m_dot_evap*solve_latent_heat_evap(T,P_tank) + m_dot_cond*solve_latent_heat_condens(T,P_tank) - P_tank*V_dot_liq + Q_dot_net
    return U_dot


def solve_T_dot(t, T, rho, m, u_prev, U_dot, rho_dot, m_dot):

    partial_dP_drho_const_T = CP.PropsSI('d(u)/d(D)|T', 'T', T, 'D', rho, 'N2O') #TESTING
    Cv = CP.PropsSI('CVMASS', 'T', T, 'D', rho, 'N2O')

    #estimate u!
    u = u_prev + (t*U_dot)/m #NOTE: I THINK THIS MAKES SENSE BUT DOUBLE CHECK

    T_dot = (1/Cv)*( (1/m) * (U_dot - u*m_dot) - (partial_dP_drho_const_T* rho_dot) )
    return T_dot



def solve_Q_dot_conduction():




def solve_thermo_properties(t, T_liq, T_gas, m_liq, m_gas, m_dot_liq, m_dot_gas, V_TANK, V_liq_prev, V_gas_prev, rho_liq_prev, rho_gas_prev):

    # Initial guess for dV/dt_liq
    V_dot_liq_guess = 0.005  # BUG CHANGE THIS --> PASS IN PREVIOUS AS INPUT You can set this based on the previous timestep or system behavior 
    V_dot_liq_iter = 0.005
    
    volume_tolerance = 1e-6 
    pressure_tolerance = 1e4
    max_iterations = 200 #TODO: define on max V_dot? can we thermodynamically figure out a maximum with some sort of hand calc assumption?

    for i in range(max_iterations):
        # Calculate dV/dt_vap based on dV/dt_liq
        V_dot_gas = -V_dot_liq_guess

        #use initial guess to solve future V_liq and V_gas with fwd euler method:
        V_liq_est = V_liq_prev + t * V_dot_liq_guess
        V_gas_est = V_gas_prev + t * V_dot_gas
        
        # Calculate d_rho/dt for liquid and vapor
        rho_dot_liq = (1 / V_liq_est) * m_dot_liq - (m_liq / V_liq_est**2) * V_dot_liq_guess 
        rho_dot_gas = (1 / V_gas_est) * m_dot_gas - (m_gas / V_gas_est**2) * V_dot_gas       
        
        # Update densities using Forward Euler method
        rho_liq = rho_liq_prev + t * rho_dot_liq
        rho_gas = rho_gas_prev + t * rho_dot_gas
        
        # Calculate pressures using Peng-Robinson EOS #NOTE: try this first, if it doesnt look good try coolprop but that might not work not sure we can see 
        pr_eos_liq = PR(T=T_liq, V=(MW/rho_liq), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
        P_liq = pr_eos_liq.P
        pr_eos_gas = PR(T=T_gas, V=(MW/rho_gas), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
        P_gas = pr_eos_gas.P
    
        
        # Check constraints:
        volume_constraint = abs(V_TANK - (m_liq / rho_liq + m_gas / rho_gas))
        pressure_constraint = abs(P_liq - P_gas)
        #TODO: pressure derivative constraint?
        
        # If both constraints are satisfied, return the solution
        if volume_constraint < volume_tolerance and pressure_constraint < pressure_tolerance:
            return P_gas, rho_liq, rho_gas, V_liq_est, V_gas_est, V_dot_liq_guess, V_dot_gas #, V_dot_liq_guess

        V_dot_liq_guess += V_dot_liq_iter
        i+=1

    raise Exception("Failed to converge")





def system_of_liq_odes(t, y, m_dot_inj):

    T_liq, T_gas, m_liq, m_gas, T_wall_liq, T_wall_gas = y  # Unpack state variables
    #NOTE: I think we are passing in previous thermo properties

    #solve thermo properties?
    P_tank, rho_liq, rho_gas, V_liq, V_gas, V_dot_liq, V_dot_gas = solve_thermo_properties(t, T_liq, T_gas, m_liq, m_gas, m_dot_liq, m_dot_gas, V_TANK, V_liq_prev, V_gas_prev, rho_liq_prev, rho_gas_prev)

    #NOTE: ORDER LIKELY BADDDD,using prev values so solving thermo properties at end?

    #solve mass transfer
    #(1) mass transfer from injector already solved
    m_dot_inj = spi_model(Cd_1, A_inj_1, P_tank, P_cc, rho_liq)
    #(2) mass transfer by condensation
    m_dot_cond = solve_m_dot_condensed(T_gas, P_tank, V_gas)
    #(3) mass transfer by evaporation
    m_dot_evap = (Q_dot_liq_to_sat_surf - Q_dot_sat_surf_to_gas) / (solve_latent_heat_evap(T_liq,P_tank)) #NOTE: P sure its T_liq but might be wrong

    m_dot_liq = -m_dot_evap + m_dot_cond - m_dot_inj
    m_dot_gas = m_dot_evap - m_dot_cond

    #NOTE: BUG BUG heat transfer signs likely wrong, I think everything should be + because Q_dot eqns will return signed vals

    #solve heat transfer terms
    T_sat = CP.PropsSI('T','P',P_tank,'Q',0,'N2O')

    # (2) from saturated surface to gas
    Q_dot_sat_surf_to_gas = solve_Q_natural_convection(T_gas, T_sat, P_tank, rho_gas, 0.15, 0.333, "N2O")
    # (3)  from liq to saturated surface
    Q_dot_liq_to_sat_surf = solve_Q_natural_convection(T_liq, T_sat, P_tank, rho_liq, 0.15, 0.333, "N2O")
    # (4) [natural convection] from liq wall to liq
    Q_dot_liq_wall_to_liq = solve_Q_natural_convection(T_liq, T_wall_liq, P_tank, rho_liq, 0.59, 0.25, "N2O")
    # (5) [natural convection] from gas wall to gas
    Q_dot_gas_wall_to_gas = solve_Q_natural_convection(T_gas, T_wall_gas, P_tank, rho_gas, 0.59, 0.25, "N2O")

    #NOTE: USE ambient properties for air, T_2 will be respective wall temperature (RK var)
    # (6) [natural convection] from atm to liq wall
    Q_dot_atm_to_liq_wall = solve_Q_natural_convection(T_ATM, T_wall_liq, P_ATM, RHO_ATM, 0.59, 0.25, "air")
    # (7) [natural convection] from atm to vap wall
    Q_dot_atm_to_gas_wall = solve_Q_natural_convection(T_ATM, T_wall_gas, P_ATM, RHO_ATM, 0.59, 0.25, "air")
    # (8) [conduction] from liq wall to vap wall
    Q_dot_liq_wall_to_gas_wall = solve_Q_dot_conduction()

    #NOTE: for T_dot_wall_liq:  IN: (6)      OUT: (4) AND (8)
    T_dot_wall_liq = ( Q_dot_atm_to_liq_wall - Q_dot_liq_wall_to_liq - Q_dot_liq_wall_to_gas_wall +  m_dot_wall*(0.15)*( T_wall_gas - T_wall_liq) ) / (0.15*m_wall_liq)

    #NOTE: for T_dot_wall_vap:  IN: (7) and (8)      OUT: (5)
    T_dot_wall_gas = ( Q_dot_atm_to_gas_wall - Q_dot_gas_wall_to_gas - Q_dot_liq_wall_to_gas_wall +  m_dot_wall*(0.15)*(T_wall_liq - T_wall_gas) ) / (0.15*m_wall_gas)
    

    #NOTE: for Q_dot_liq:   IN:(2)       OUT: (3)
    U_dot_liq = solve_U_dot(T_liq, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_liq , (Q_dot_liq_wall_to_liq - Q_dot_liq_to_sat_surf) )

    #NOTE: for Q_dot_liq:   IN:(5)       OUT: (NONE)
    U_dot_gas = solve_U_dot(T_gas, P_tank, m_dot_inj, m_dot_evap, m_dot_cond, V_dot_gas, Q_dot_gas_wall_to_gas)


    T_dot_liq = solve_T_dot(t, T_liq, rho_liq, m_liq, u_prev_liq, U_dot_liq, rho_dot_liq, m_dot_liq)

    T_dot_gas = solve_T_dot(t, T_gas, rho_gas, m_gas, u_prev_gas, U_dot_gas, rho_dot_gas, m_dot_gas)


    
    return [T_dot_liq, T_dot_gas, m_dot_liq, m_dot_gas, T_dot_wall_liq, T_dot_wall_gas]








class model():

    #NOTE/TODO: THIS IS PROBABLY THE SAME as bens, change if not
    def __init__(self, TIMESTEP, T_atm, m_nos, Cd_1, A_inj_1, V_tank, Diam_tank, P_tank, P_cc, all_error, inj_model):
        
        #injector input
        self.Cd_1 = Cd_1
        self.P_cc = P_cc
        self.A_inj_1 = A_inj_1

        #NOTE: Setup - Assuming Tank Starts at Sat Conditions
        self.TIMESTEP = TIMESTEP
        self.convergence_percent = 0.001 #NOTE: ADJUST AS REQUIRED,

        self.V_tank = V_tank
        self.m_nos = m_nos
        v_init = V_tank/m_nos  # Initial volume in m^3

        rho_sat_vap = CP.PropsSI('D', 'P', P_tank, 'Q', 1, "N2O")
        rho_sat_liq = CP.PropsSI('D', 'P', P_tank, 'Q', 0, "N2O")
        rho_sat_tank = 1/v_init

        x_tank = ( (1/rho_sat_tank)-(1/rho_sat_liq) ) / ( (1/rho_sat_vap)-(1/rho_sat_liq) )

        ### initial cv values from initial sat liq vap state

        #vapor cv
        self.m_vap = x_tank*self.m_nos

        self.v_vap = 1/rho_sat_vap
        self.V_vap = self.v_vap*self.m_vap
        self.P_vap = P_tank #+ 100 #NOTE: TESTING PERTURB

        self.P_dot_vap = 0
        self.rho_dot_vap = 0

        #liquid cv
        self.m_liq = self.m_nos-self.m_vap

        self.rho_liq = rho_sat_liq
        self.v_liq = 1/self.rho_liq
        self.V_liq = self.V_tank-self.V_vap
        self.P_liq = P_tank + 100 #NOTE: TESTING PERTURB

        """
        pr_eos = PR(P=self.P_vap, V= (MW*self.v_vap), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
        self.T_vap = pr_eos.solve_T(self.P_vap, (MW*self.v_vap) )
        """

        self.T_vap = CP.PropsSI('T', 'P', P_tank, 'Q', 1, "N2O")  #NOTE: TESTING PERTURB


        self.rho_dot_liq = 0
        self.P_dot_liq = 0
        self.V_dot_liq = 0
        self.T_liq = self.T_vap #starting at equillibrium

        #initially at equillibrium
        self.Q_dot_evap = 0
        self.m_dot_evap = 0
        self.m_dot_inj = 0 
                
        self.P_vap_prev = P_tank
        self.P_liq_prev = P_tank

        self.rho_vap_prev = (1/self.v_vap)
        self.rho_liq_prev = (1/self.v_liq)

        self.rho_exit = 1/self.v_liq

        self.t = 0

        print(self.V_liq, self.T_liq, self.rho_liq, self.P_liq)
        

        

        



    def inst(self, P_cc):

        self.t+=self.TIMESTEP
        
        if self.m_liq >= 0: #liquid phase!
            
            #enter 6 ode system




            return 1
            


        else: #vapor phase

            return 0










t = 0
TIMESTEP = 1e-4
T_atm = None #not used rn
m_nos = 20
Cd_1 = 0.425
A_inj_1 = 0.00003 #m^3 NOTE: GUESS
V_tank = 0.0354
Diam_tank = None #NOTE: fine for testing but we will need to pass this in in a bit
P_tank = 45e5
P_cc = 1.03e6
all_error = None #unused
inj_model = None #TODO: implement


tank = model(TIMESTEP, T_atm, m_nos, Cd_1, A_inj_1, V_tank, Diam_tank, P_tank, P_cc, all_error, inj_model)
while(t<5):
    tank.inst(P_cc)
    t+=TIMESTEP
    print(t, tank.P_liq, tank.V_liq,"\n")
    #print("\n")