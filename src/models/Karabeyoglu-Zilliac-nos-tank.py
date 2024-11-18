
from scipy.optimize import minimize
from thermo import Chemical
from thermo import PR
import numpy as np
from scipy.integrate import solve_ivp
import rocketprops
import CoolProp.CoolProp as CP

# Global Constants:
R_U = 8.31446 #J/(mol K) 

#T_REF = 298.15 #K
#P_REF = 101325 #Pa
#NOTE: do i even use this? ^
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






# Heat Transfer Functions and LHS

def latent_heat_vap(T): 
    return (CP.PropsSI('H', 'T', T, 'Q', 1, "N2O") - CP.PropsSI('H', 'T', T, 'Q', 0, "N2O")) #J/kg
#BUG: assuming equilibrium here, that seems wrong... #NOTE: YES BAD !!!

def solve_Q_dot_evap(T, m_dot_evap):
    #print("m_dot_evap", m_dot_evap)
    return m_dot_evap*latent_heat_vap(T) #J/s

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
    m_dot_evaporated_liq = h_lsb * (CS_AREA/latent_heat_vap(T_1_l))*(T_1_l-T_1_s)
    return m_dot_evaporated_liq


def calculate_LHS(P_1_l, T_1_l, rho_1_l, m_dot_inj):

    h_l = CP.PropsSI('H', 'T', T_1_l, 'P', P_1_l, 'N2O')

    m_dot_evap = solve_m_dot_evaporated_liq(T_1_l, P_1_l, rho_1_l)
    Q_dot_evap = solve_Q_dot_evap(T_1_l, m_dot_evap)

    h_evap = latent_heat_vap(T_1_l)

    Q_dot_ls = solve_Q_dot_ls(T_1_l, P_1_l, rho_1_l)


    #print(f"LHS: Q_dot_evap={Q_dot_evap} (J/s), Q_dot_ls={Q_dot_ls} (J/s), m_dot={m_dot} (kg/s), h_l={h_l} (J/kg), m_dot_evap={m_dot_evap}, h_evap={h_evap} (J/kg)\n")
    #print("checking LHS SIGN", (Q_dot_evap -Q_dot_ls), (m_dot*h_l), ( m_dot_evap*h_evap) )
    
    lhs = (Q_dot_evap -Q_dot_ls) + m_dot_inj*h_l + m_dot_evap*h_evap
    return lhs








# Thermodynamic Functions and RHS

def solve_A_B_Z_alpha(T, rho):
    #convert rho to V_m
    V_m = MW / rho

    pr_eos_vap = PR(T=T, V=V_m, Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega) #TODO: check units
    P = pr_eos_vap.P #allegedly this is how pressure is solved (on init?)

    Z = (P*V_m)/(R_U*T)

    A = (pr_eos_vap.a*P)/( ((R_U/MW)**2)*T**2)
    B = (pr_eos_vap.b*P)/((R_U/MW)*T)

    alpha = (1+ KAPPA * (1 - np.sqrt(T/n2o.T_cr))) 

    return A, B, Z, alpha

def solve_F_1(T, rho):


    A, B, Z, alpha = solve_A_B_Z_alpha(T, rho)

    F_1 = ((R_U*T)/MW) * (np.log( (Z+2.414*B)/(Z-0.414*B) ) * (A/(5.657*B)) * (KAPPA/(n2o.Tc*alpha)) * (np.sqrt(alpha/T_REF) + KAPPA) )
    return F_1


def solve_F_2(T, rho):
    A, B, Z, alpha = solve_A_B_Z_alpha(T, rho)

    F_2 = ((R_U*T)/MW) * ( (-1)*(Z/rho) * (A/((Z**2)+2*B*Z-(B**2))) * (1+KAPPA*np.sqrt(T_REF/alpha)) )
    return F_2


def calculate_RHS(P:float, T:float, rho:float, m:float, P_dot:float, T_dot:float, rho_dot:float, m_dot:float, V_dot:float):

    F_1 = solve_F_1(T, rho)
    F_2 = solve_F_2(T, rho)

    cv_ig = (n2o.HeatCapacityGas.T_dependent_property(T)/MW) - (R_U/MW)

    u = CP.PropsSI('U', 'T', T, 'P', P, 'N2O') #J/kg


    RHS = (P*V_dot + P_dot*(m/rho)) + m_dot*u + m*(F_1 + cv_ig)*T_dot + m*F_2*rho_dot
    return RHS









### V_dot eqns for explicit property estimation (recall vapor phase V_dot = 0)

def liq_phase_v_dot_liq(V, rho, rho_dot, m_dot_inj, m_dot_evap):
    v_dot = ( -m_dot_inj -m_dot_evap -V*rho_dot) / rho
    return v_dot

#do i even need this? just sign result of other?
def liq_phase_v_dot_vap(V_liq, rho_liq, rho_dot_liq, m_dot_inj, m_dot_evap):
    return -liq_phase_v_dot_liq(V_liq, rho_liq, rho_dot_liq, m_dot_inj, m_dot_evap)




### T_dot eqns for explicit property estimation

def vap_phase_T_dot(T, rho, rho_dot, m, m_dot_prop, n2o, ke, Q_dot_gas_wall):

    #actually likely makes sense to keep T separate or define a new PR object????
    #BUG: this is definitely wrong, want T_sat?
    pr_eos_vap = PR(T=T, V=(MW/rho), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega) #TODO: check units
    P = pr_eos_vap.P #allegedly this is how pressure is solved (on init?)

    F_1 = solve_F_1(P, T, rho, n2o)
    F_2 = solve_F_2(P, T, rho, n2o)


    ####


    cv_ideal_gas = n2o.Cp_ideal_gas_mass()  - (R_U/MW)

    T_dot = (Q_dot_gas_wall - m_dot_prop*( (P/rho) + 0.5*(ke**2) ) - m*F_2*rho_dot) / (m*(F_1+cv_ideal_gas))
    return T_dot


def liq_phase_v_dot_liq(P, T, rho, Q_dot_net, m_dot_inj, m_dot_evap, h_evap, m):

    F_1 = solve_F_1(T, rho)
    F_2 = solve_F_2(T, rho)

    u = CP.PropsSI('U', 'T', T, 'P', P, 'N2O') #J/kg
    h = CP.PropsSI('U', 'T', T, 'P', P, 'N2O') #J/kg

    cv_ideal_gas = n2o.Cp_ideal_gas_mass() - (R_U/MW)

    T_dot = (Q_dot_net - m_dot_inj*h - m_dot_evap*h_evap -(P*V_dot + P_dot*V) - m_dot_inj*u - m*F_2*rho_dot ) / (m*(F_1+cv_ideal_gas))


    return T_dot


def liq_phase_T_dot_vap(P,T,rho,):

    F_1 = solve_F_1(P, T, rho, n2o)
    F_2 = solve_F_2(P, T, rho, n2o)

    u = CP.PropsSI('U', 'T', T, 'P', P, 'N2O') #J/kg
    h = CP.PropsSI('U', 'T', T, 'P', P, 'N2O') #J/kg

    cv_ideal_gas = n2o.Cp_ideal_gas_mass() - (R_U/MW)

    T_dot = (Q_dot_net - m_dot_evap*h_evap -(P*V_dot + P_dot*V) - m_dot_inj*u - m*F_2*rho_dot ) / (m*(F_1+cv_ideal_gas))


    return T_dot







class model():

    #NOTE/TODO: THIS IS PROBABLY THE SAME as bens, change if not
    def __init__(self, oxidizer, TIMESTEP, T_atm, m_nos, Cd_1, A_inj_1, V_tank, Diam_tank, P_tank, P_cc, all_error, inj_model):
        
        ###Injector Constants


        # setup - start by assuming the tank is in thermal equillibrium
        
        rho_sat_vap = CP.PropsSI('D', 'P', P_tank, 'Q', 1, "N2O")
        rho_sat_liq = CP.PropsSI('D', 'P', P_tank, 'Q', 0, "N2O")
        rho_sat_tank = m_nos/V_tank

        x_tank = ( (1/rho_sat_tank)-(1/rho_sat_liq) ) / ( (1/rho_sat_vap)-(1/rho_sat_liq) )

        self.m_vap = x_tank*m_nos
        self.m_liq = m_nos-self.m_vap

        self.v_vap = 1/rho_sat_vap
        self.P_vap = P_tank
        self.T_vap = T_atm
    

        self.v_liq = 1/rho_sat_liq
        self.P_liq = P_tank
        self.T_liq = T_atm


        self.Q_dot_evap = 0
        self.m_dot_evap = 0
        self.m_dot_inj = 0
        self.h_evap = (CP.PropsSI('H', 'T', self.T_vap, 'Q', 1, "N2O") - CP.PropsSI('H', 'T', self.T_liq, 'Q', 0, "N2O")) #J/kg
        
        self.P_vap_prev = P_tank
        self.P_liq_prev = P_tank
        
        self.convergence_percent = 0.001 #NOTE: ADJUST AS REQUIRED, better yet make it a percent
        

        



    def inst(self, P_cc):

        self.m_dot_inj = spi_model(self.Cd_1, self.A_inj_1, self.P_vap, P_cc, self.rho_exit)
        #NOTE: this should be vapor or liquid pressure? solve in init?
        
        if self.m_liq >= 0: #two phases in the tank!

            #recall m_dot, Q_dot, h_evap terms stored as attributes from

            sol = solve_ivp(system_of_thermo_property_odes, TIMESTEP, y0, method='RK45')
            






            RHS = calculate_RHS(self.P_liq, self.T_liq, (1/self.v_liq), self.m_liq, 0, 0, 0, self.m_dot_inj)
            LHS = calculate_LHS(self.P_liq, self.T_liq, (1/self.v_liq), self.m_dot_inj)
            convergence_criteria = (np.abs(RHS) + np.abs(LHS) ) / 2 * self.convergence_percent

            while( convergence_criteria < ((RHS-LHS)**2) ):

                #solve mass transfer for both cv

                #solve m_dot evap
                #NOTE: i think this is the wrong formula, this will depend on both CV
                m_dot_evap = solve_m_dot_evaporated_liq(self.T_liq, self.P_liq, (1/self.v_liq))          


                #solve heat transfer for both cv

                Q_dot_evap = solve_Q_dot_evap(self.T_dot_liq, m_dot_evap)
                Q_dot_ls = solve_Q_dot_ls(self.T_liq, self.P_liq, (1/self.v_liq)) 


                #solve T_dot and rho_dot #NOTE: missing a lot of steps here

                #this is most defnitely a runge kutta method!


                ### solve volume and temperature change



                #PR EOS to solve new pressure for both cv NOTE: takes molar volume as input!!!!
                pr_eos_liq = PR(T=self.T_liq, V=(MW*self.v_liq),  Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
                self.P_liq = pr_eos_liq.P
                P_dot_liq = (self.P_liq - self.P_liq_prev) / self.TIMESTEP

                pr_eos_vap = PR(T=self.T_vap, V=(MW*self.v_vap),  Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
                self.P_vap = pr_eos_vap.P
                P_dot_vap = (self.P_vap - self.P_vap_prev) / self.TIMESTEP
        




                #solved new T,P,rho, for both cv check if density converged in convergence loop
                #resolve LHS, RHS
                RHS = calculate_RHS(self.P_liq, self.T_liq, (1/self.v_liq), self.m_liq, P_dot_liq, T_dot_liq, self.v_dot_liq, (self.m_dot_inj+m_dot_evap) )
                LHS = calculate_LHS(self.P_liq, self.T_liq, (1/self.v_liq), self.m_dot_inj)
            #NOTE: I DONT THINK THIS WOULD ALLOW VAPOR CV TO CONTROL CONVERGENCE
                













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