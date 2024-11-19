from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical

# Global Constants:
T_REF = 298.15 #K
P_REF = 101325 #Pa

R_U = 8.31446 #J/(mol K) 

n2o = Chemical('nitrous oxide')

MW = (n2o.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o.omega - 0.26992*n2o.omega**2

TANK_DIAM = 0.0254*5.5 #m
CS_AREA = 0.25*np.pi*(TANK_DIAM**2) #m^2
g = 9.81 #m/s^2



#TODO: update with other injector model once we get this thing up
def spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit):
    m_dot_spi = Cd_hem_spi_dyer * A_inj_ox * np.sqrt( 2 * rho_tank_exit * (P_1 - P_2)  )
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



def solve_A_B_Z_alpha(T, rho):
    #convert rho to V_m
    V_m = MW / rho

    pr_eos_vap = PR(T=T, V=V_m, Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega) #TODO: check units
    P = pr_eos_vap.P #allegedly this is how pressure is solved (on init?)

    Z = (P*V_m)/(R_U*T)

    A = (pr_eos_vap.a*P)/( ((R_U/MW)**2)*T**2)
    B = (pr_eos_vap.b*P)/((R_U/MW)*T)

    alpha = (1+ KAPPA * (1 - np.sqrt(T/n2o.Tc))) 

    return A, B, Z, alpha

def solve_F_1(T, rho):
    A, B, Z, alpha = solve_A_B_Z_alpha(T, rho)

    F_1 = ((R_U*T)/MW) * (np.log( (Z+2.414*B)/(Z-0.414*B) ) * (A/(5.657*B)) * (KAPPA/(n2o.Tc*alpha)) * (np.sqrt(alpha/T_REF) + KAPPA) )
    return F_1


def solve_F_2(T, rho):
    A, B, Z, alpha = solve_A_B_Z_alpha(T, rho)

    F_2 = ((R_U*T)/MW) * ( (-1)*(Z/rho) * (A/((Z**2)+2*B*Z-(B**2))) * (1+KAPPA*np.sqrt(T_REF/alpha)) )
    return F_2






def solve_liq_phase_V_dot_liq(V, rho, rho_dot_liq, m_dot_inj, m_dot_evap):

    V_dot_liq = (-m_dot_inj - m_dot_evap - V*rho_dot_liq) / rho #NOTE: NOT SPECIFIC VOLUME!!!!
    print("V_dot_liq", V_dot_liq)
    return V_dot_liq


    
def solve_liq_phase_T_dot_liq(P, T, rho, Q_dot_net, m_dot_inj, m_dot_evap, h_evap, m, P_dot_liq, V_dot_liq, rho_dot_liq):

    F_1 = solve_F_1(T, rho)
    F_2 = solve_F_2(T, rho)

    u = CP.PropsSI('U', 'T', T, 'P', rho, 'N2O') #J/kg
    h = CP.PropsSI('U', 'T', T, 'P', P, 'N2O') #J/kg

    cv_ideal_gas = (n2o.HeatCapacityGas.T_dependent_property(T)/MW) - (R_U/MW)

    #NOTE: FIRST TIMESTEP STARTS IN EQUILIBRIUM
    V_dot_liq = 0
    P_dot_liq = 0
    rho_dot_liq = 0

    V = m/rho

    T_dot = (Q_dot_net - m_dot_inj*h - m_dot_evap*h_evap -(P*V_dot_liq + P_dot_liq*V) - m_dot_inj*u - m*F_2*rho_dot_liq ) / (m*(F_1+cv_ideal_gas))

    print("T_dot", T_dot)
    return T_dot

def solve_liq_phase_T_dot_vap(t, y, P, rho, Q_dot_net, m_dot_inj, m_dot_evap, h_evap, m, P_dot_vap, V_dot_vap, rho_dot_vap):

    T = y

    F_1 = solve_F_1(T, rho)
    F_2 = solve_F_2(T, rho)

    u = CP.PropsSI('U', 'T', T, 'P', rho, 'N2O') #J/kg

    cv_ideal_gas = (n2o.HeatCapacityGas.T_dependent_property(T)/MW) - (R_U/MW)

    V = m/rho

    m_dot_vap = m_dot_evap #I am neglecting condensation for now, no tank pressurant gas!
    T_dot = (Q_dot_net + m_dot_evap*h_evap -(P*V_dot_vap + P_dot_vap*V) - m_dot_inj*u - m*F_2*rho_dot_vap ) / (m*(F_1+cv_ideal_gas)) #NOTE: TAKE CARE W SIGNS HERE

    return T_dot






#TODO: what is the best way to handle the change? check states with coolprop??????
def solve_latent_heat_evap(T,P): 
    
    phase = CP.PropsSI('Phase', "T", T, "P", P, "N2O")

    print("phase: ", phase)

    h_delta = (CP.PropsSI('H', 'T', T, 'Q', 1, "N2O") - CP.PropsSI('H', 'T', T, 'Q', 0, "N2O")) #J/kg

    if phase == 0: #subcooled fluid 
        h_evap = h_delta + (CP.PropsSI('H', 'T', T, 'Q', 0, "N2O") - CP.PropsSI('H', 'T', T, 'P', P, "N2O") ) #J/kg

    if phase == 6: #sat liq vapor
        h_evap = h_delta

    return h_evap #J/kg


def solve_Q_dot_evap(T, P, m_dot_evap):
    #print("m_dot_evap", m_dot_evap)
    return m_dot_evap*solve_latent_heat_evap(T, P) #J/s

def solve_h_ls(T_1_l, P_1_l, rho_1_l, T_1_s):
    #NOTE: THIS MIGHT NOT APPLY TO NOS
    C_ls = 0.27
    n_ls = 0.25

    k_l = get_thermal_conductivity(T_1_l, P_1_l) #(W/(m K))
    Cp_l = CP.PropsSI('C', 'T', T_1_l, 'P', P_1_l, 'N2O') #J/(kg K)
    visc_l = get_viscosity(T_1_l, P_1_l) #(Pa s)

    dV_dT_P = n2o.VolumeLiquid.TP_dependent_property_derivative_T(T_1_l, P_1_l)  #TODO: CHECK THIS
    #print(f"\n dV_dT_P = {dV_dT_P}")
    beta = dV_dT_P*rho_1_l 

    Gr = ((TANK_DIAM**3)*(rho_1_l**2)*g*beta*np.abs(T_1_l - T_1_s) ) / (visc_l**2)
    Pr = (Cp_l*visc_l)/ k_l

    X_ls = Gr*Pr
    h_ls = C_ls * (k_l/TANK_DIAM) * X_ls**n_ls

    return h_ls

#Q_dot_LS - convection between bulk liquid surface and liquid surface layer
def solve_Q_dot_ls(T_1_l, P_1_l, rho_1_l):

    T_1_s = CP.PropsSI('T', 'P', P_1_l, 'Q', 0, 'N2O')
    h_ls = solve_h_ls(T_1_l, P_1_l, rho_1_l, T_1_s)

    Q_dot_ls = h_ls*CS_AREA*(T_1_l - T_1_s) #liquid temp minus surface temp --> assuming surface temp is sat temp?

    #print("\nCHECKING Q_DOT_LS: ", Q_dot_ls, T_1_l, T_1_s, h_ls, "\n")
    return Q_dot_ls

def solve_m_dot_evaporated_liq(T_1_l, P_1_l, rho_1_l):
    #NOTE: fluid restricted to one component and surface temperature assumed to equal saturation temperature
    
    T_1_s = CP.PropsSI('T', 'P', P_1_l, 'Q', 0, 'N2O')
    h_ls = solve_h_ls(T_1_l, P_1_l, rho_1_l, T_1_s)
    
    E = 2.1e4 #correction coeff for N2O from thesis
    h_lsb = E*h_ls
    
    #latent heat from ideal gas relation under model? #     T_1 and T_Liq same?
    m_dot_evaporated_liq = h_lsb * (CS_AREA/solve_latent_heat_evap(T_1_l, P_1_l))*(T_1_l-T_1_s)
    return m_dot_evaporated_liq







def system_of_liq_odes(t, y, P_tank, m_dot_inj, m_dot_evap, Q_dot_net, m, rho, rho_dot_liq):

    T, v = y  # Unpack state variables

    #print("before the crash: ", T, P_tank, v)


    h_evap = (CP.PropsSI('H', 'T', T, 'D', rho, "N2O") - CP.PropsSI('H', 'T', T, 'D', rho, "N2O")) #J/kg

    V_dot_liq = solve_liq_phase_V_dot_liq( (v*m), rho, rho_dot_liq, m_dot_inj, m_dot_evap) 

    P_dot_liq = 0 #TODO: PASS THIS IN from previous P_dot_liq

    print("FIRST TIMESTEP T_DOT_LIQ: ", P_tank, T, rho, Q_dot_net, m_dot_inj, m_dot_evap, h_evap, m, P_dot_liq, 0, rho_dot_liq)
    T_dot_liq = solve_liq_phase_T_dot_liq(P_tank, T, rho, Q_dot_net, m_dot_inj, m_dot_evap, h_evap, m, P_dot_liq, 0, rho_dot_liq) #NOTE: just to use previous val V_dot_liq = 0
    
    print(V_dot_liq, T_dot_liq)
    return [T_dot_liq, V_dot_liq]








class model():

    #NOTE/TODO: THIS IS PROBABLY THE SAME as bens, change if not
    def __init__(self, oxidizer, TIMESTEP, T_atm, m_nos, Cd_1, A_inj_1, V_tank, Diam_tank, P_tank, P_cc, all_error, inj_model):
        
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
        self.P_vap = P_tank

        self.P_dot_vap = 0
        self.rho_dot_vap = 0

        #liquid cv
        self.m_liq = self.m_nos-self.m_vap

        self.v_liq = 1/rho_sat_liq
        self.V_tank = self.V_tank-self.V_vap
        self.P_liq = P_tank

        pr_eos = PR(P=self.P_vap, V= (MW*self.v_vap), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
        self.T_vap = pr_eos.solve_T(self.P_vap, (self.MW*self.v_vap) )

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

        self.rho_vap_prev = 1/self.v_vap
        self.rho_liq_prev = 1/self.v_liq

        
        

        

        



    def inst(self, P_cc):

        self.m_dot_inj = spi_model(self.Cd_1, self.A_inj_1, self.P_vap, P_cc, self.rho_exit)
        #NOTE: this should be vapor or liquid pressure? solve in init?
        
        if self.m_liq >= 0: #liquid phase!

            # Initial state setup
            y0_liq = [self.T_liq, (self.v_liq*self.m_liq)]

            # Solve the system

            self.m_liq-=self.m_dot_inj*self.TIMESTEP

            Q_dot_net_liq = self.Q_dot_evap #solve simple case to start
            solution = solve_ivp(system_of_liq_odes, [0, self.TIMESTEP], y0_liq, args=(self.P_liq, self.m_dot_inj, self.m_dot_evap, Q_dot_net_liq, self.m_liq, (1/self.v_liq), self.rho_dot_liq),
                                method='RK45', rtol=1e-8, atol=1e-8)

            # Extract results
            y_liq = solution.y

            self.T_liq = y_liq[0, -1]
            self.V_liq = y_liq[1, -1]

            pr_eos = PR(T=self.T_liq, V= (MW*self.V_liq/self.m_liq), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
            self.P_liq = pr_eos.P




            #### NOW RESOLVE LIQ TIME DERIVATIVES FOR VAPOR PHASE AND NEXT RUNGE KUTTA STEP


            self.rho_liq = self.m_liq/self.V_liq

            self.m_dot_evap = solve_m_dot_evaporated_liq(self.T_liq, self.P_liq, self.rho_liq)
            self.Q_dot_evap = solve_Q_dot_evap(self.T_liq, self.P_liq, self.m_dot_evap)

            #NOTE: P_vap? P_liq?
            self.m_dot_inj = spi_model(self.Cd_1, self.A_inj_1, self.P_vap, self.P_cc, self.rho_liq) #TODO: FILL IN FOR TEST CASE THEY USED SPI AS WELL

            self.rho_dot_liq = (self.rho_liq - self.rho_liq_prev)/self.TIMESTEP
            self.V_dot_liq = solve_liq_phase_V_dot_liq(self.V_liq, self.rho_liq, self.rho_dot_liq, self.m_dot_inj, self.m_dot_evap)

            self.h_evap = solve_latent_heat_evap(self.T_liq, self.P_liq)

            Q_dot_net_liq = self.Q_dot_evap #simplified model just to get it to work at the start


            #print("INTO NEXT TIMESTEP T_DOT_LIQ: ", P_liq, T_liq, rho_liq, Q_dot_net_liq, m_dot_inj, m_dot_evap, h_evap, m_liq, P_dot_liq, V_dot_liq, rho_dot_liq)
            self.T_dot_liq = solve_liq_phase_T_dot_liq(self.P_liq, self.T_liq, self.rho_liq, Q_dot_net_liq, self.m_dot_inj, self.m_dot_evap, self.h_evap, self.m_liq, self.P_dot_liq, self.V_dot_liq, self.rho_dot_liq)
            self.P_dot_liq = (self.P_liq - self.P_liq_prev)/self.TIMESTEP


            ### VAPOR CV
            self.V_dot_vap = (-self.V_dot_liq) #--> use this to solve rho_vap below!!!!!

            #fwd euler step to solve rho_vap
            #print("checking sign convention: ", self.v_vap)
            self.V_vap += self.V_dot_vap*self.TIMESTEP
            #print("checking sign convention after: ", self.v_vap)

            Q_dot_net_vap = self.Q_dot_evap #simplified model just to get it to work at the start

            solution_vap = solve_ivp(solve_liq_phase_T_dot_vap, [0, self.TIMESTEP], [self.T_vap], args=(self.P_vap, 1/(self.v_vap), Q_dot_net_vap, self.m_dot_inj, self.m_dot_evap, self.h_evap, self.m_vap, self.P_dot_vap, self.V_dot_vap, self.rho_dot_vap),
                                method='RK45', rtol=1e-8, atol=1e-8)

            # Extract results
            y_vap = solution_vap.y
            self.T_vap = y_vap[0, -1]
            self.V_vap = self.V_tank - self.V_liq

            pr_eos = PR(T=self.T_vap, V= (MW*self.V_vap/self.m_vap), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
            self.P_vap = pr_eos.P

            ###resolve derivatives
            self.P_dot_vap = (self.P_vap - self.P_vap_prev)/self.TIMESTEP
            #NOTE: already have V_dot_vap above 
            self.rho_vap = self.m_vap/self.V_vap
            self.rho_dot_vap = (self.rho_vap - self.rho_vap_prev)/self.TIMESTEP

            self.T_dot_vap = solve_liq_phase_T_dot_vap(None, self.T_vap, self.P_vap, 1/(self.v_vap), Q_dot_net_vap, self.m_dot_inj, self.m_dot_evap, self.h_evap, self.m_vap, self.P_dot_vap, self.V_dot_vap, self.rho_dot_vap),





            #print("after", self.T_vap, self.P_vap, self.rho_vap, self.V_vap)

            #print("vap derivatives: ", self.T_dot_vap, self.V_dot_vap, self.rho_dot_vap, self.P_dot_vap)






        #else: #vapor phase
            














        else: #only gas in tank (only one cv same idea though?)


























"""#OLD:
            #update conservation of mass

            self.m_vap -= self.TIMESTEP * self.m_dot

            rho_vap_prev = self.rho_vap
            self.rho_vap = self.m_vap/self.V_tank

            rho_dot = (self.rho_vap - rho_vap_prev) / self.TIMESTEP

            #solve Q_dot_vap_wall
            Q_dot_vap_wall = solve_Q_dot_vap_wall()
            
            #setting ke term to 0 for right now, likely negligible
            u_e = 0

            #NOTE: solve_ivp returns a object with a bunch of info in it
            sol = solve_ivp( lambda T: vap_phase_T_dot(T, self.rho_vap, rho_dot, self.m_vap, self.m_dot, self.n2o, u_e, Q_dot_vap_wall)
                            , self.TIMESTEP, self.T_vap, method='RK45', t_eval=self.TIMESTEP)

            self.T_vap = sol.y[0,-1]

            #now that we solved T and have rho, can update tank pressure with Peng Robinson EOS
        
            #convert rho to V_m
            V_m = (self.n2o.MW/1000) / self.rho_vap #note converting n2o.MW g/mol to kg/mol

            # Initialize the PR EOS with the given temperature and molar volume
            self.pr_eos_vap = PR(T=self.T_vap, V=V_m, Tc=self.n2o.Tc, Pc=self.n2o.Pc, omega=self.n2o.omega)
            self.P_vap = self.pr_eos_vap.P

            #resolve all thermodynamic properties

            self.rho_exit = self.rho_vap
"""