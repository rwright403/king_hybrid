from thermo import Chemical
from thermo.eos import PR

# Create a Fluid object for nitrous oxide
nitrous_oxide = Chemical('nitrous oxide')

from scipy.integrate import solve_ivp
import numpy as np

#global constants:
R_u = 8.31446 #J/K/mol

T_ref = 298.15 #K
P_ref = 101325 #Pa


def spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit):
    m_dot_spi = Cd_hem_spi_dyer * A_inj_ox * np.sqrt( 2 * rho_tank_exit * (P_1 - P_2)  )
    #print(rho_tank_exit)
    return m_dot_spi

### property relations
    
def dynamic_visc_vap(T):
    A = 2.1150e-6
    B = 0.46420
    C = 305.70

    T_max = 1000 #K
    T_min = 182 #K

    if (T < T_min) or (T > T_max):
        raise ValueError(f"Input T={T} is out of bounds! dynamic_visc_vap only valid for [{T_min},{T_max}]")
    else:
        val = A+B*T+C*T**2
        return val #kg/(m s)

def beta_liq(T):
    A = 2.781
    B = 0.27244
    C = 309.57
    D = 0.2882

    T_max = 310 #K
    T_min = 160 #K

    if (T < T_min) or (T > T_max):
        raise ValueError(f"Input T={T} is out of bounds! beta_liq only valid for [{T_min},{T_max}]")
    else:
        val = (-D/C)*np.ln(B)*(1-T/C)**(D-1)
        return val #1/K

def T_sat_vap(P):
    A = 4.80716087
    B = 967.819748
    C = 19.6368887

    val = B/(A-np.log10(P/100000))-C
    return val #K

def P_sat_vap(T):
    A = 4.80716087
    B = 967.819748
    C = 19.6368887

    T_max = 310 #K
    T_min = 140 #K

    if (T < T_min) or (T > T_max):
        raise ValueError(f"Input T={T} is out of bounds! P_sat_vap only valid for [{T_min},{T_max}]")
    else:
        val = 100000*(10**(A-(B/(T+C))))
        return val #Pa TODO: CHECK UNITS I THINK ITS IN KPA

def dynamic_visc_liq(T):
    A = 0.001877085
    B = 1.1864e-5
    C = 1.928e-8

    T_max = 300 #K
    T_min = 230 #K

    if (T < T_min) or (T > T_max):
        raise ValueError(f"Input T={T} is out of bounds! dynamic_visc_liq only valid for [{T_min},{T_max}]")
    else:
        val = A+B*T+C*T**2
        return val #kg/(m s)

def k_liq(T):
    A = 0.39815
    B = -0.0011510

    T_max = 310 #K
    T_min = 182 #K

    if (T < T_min) or (T > T_max):
        raise ValueError(f"Input T={T} is out of bounds! k_liq only valid for [{T_min},{T_max}]")
    else:
        val = A+B*T
        return val #W/(m K)

def k_vap(T):
    A = -0.007037
    B = 0.0000823

    T_max = 430 #K
    T_min = 190 #K

    if (T < T_min) or (T > T_max):
        raise ValueError(f"Input T={T} is out of bounds! k_vap only valid for [{T_min},{T_max}]")
    else:
        val = A+B*T
        return val #W/(m K)

def cp_sat_vap(T):
    A = -5956.82087
    B = 59029.4538
    C = -215342.983
    D = 276450.549
    E = 23.5743297

    T_max = 310 #K
    T_min = 150 #K

    if (T < T_min) or (T > T_max):
        raise ValueError(f"Input T={T} is out of bounds! cp_sat_vap only valid for [{T_min},{T_max}]")
    else:
        val = A + B*T + C*(T**2) + D*(T**3) + E*(T**4)
        return val #J/(kg K)

def cp_ideal_gas(T):
    A = 21.62
    B = 78.81
    C = -57.78
    D = 18.3

    T_max = 310 #K
    T_min = 150 #K

    if (T < T_min) or (T > T_max):
        raise ValueError(f"Input T={T} is out of bounds! cp_ideal_gas only valid for [{T_min},{T_max}]")
    else:
        val = A + B*T + C*(T**2) + D*(T**3) + E*(T**4)
        return val #J/(kg K)
    
def cp_sat_liq(T):
    A = -131.55551
    B = 3013.03128
    C = -14290.1471
    D = 22239.8432

    T_max = 390 #K
    T_min = 183 #K

    if (T < T_min) or (T > T_max):
        raise ValueError(f"Input T={T} is out of bounds! cp_sat_liq only valid for [{T_min},{T_max}]")
    else:
        val = A + B*T + C*(T**2) + D*(T**3) + E*(T**4)
        return val #J/(kg K)

def latent_heat_liq_vap(T,MW):
    A = 2.686e7
    B = 0.182
    C = 0.9387
    D = -0.706

    T_max = 390 #K
    T_min = 183 #K

    #NOTE: thesis includes T_crit = 309.6 K here, not sure the significance rn

    if (T < T_min) or (T > T_max):
        raise ValueError(f"Input T={T} is out of bounds! latent_heat_liq_vap only valid for [{T_min},{T_max}]")
    else:
        val = (A*(1-T_ref)**(B+C*T_ref+D*(T_ref**2)+E*(T_ref**3)))/MW 
        return val #J/(kg K) 



#Q_dot_GW - convection
def solve_Q_dot_vap_wall(T_tank_gas,T_wall_gas,rho_gas,L_scale,diam_tank,C_gas_wall,n_gas_wall):
    
    A_gas_wall = np.pi*diam_tank*L_scale #m^2

    T_gas_wall_film = (T_tank_gas + T_wall_gas)/2 #K 

    beta_gas_wall_film = beta_ ??? (T_gas_wall_film)
    dynamic_visc_vap = dynamic_visc_vap(T_gas_wall_film)
    cp_vap = ??? (T_gas_wall_film)
    k_gas_wall_film = k_vap(T_gas_wall_film)


    X_gas_wall = ((L_scale**3)*(rho_gas**2)*alpha*beta_gas_wall_film* np.abs(T_tank_gas - T_wall_gas)/(dynamic_visc_vap**2) ) * (cp_vap * dynamic_visc_vap/k_gas_wall_film)

    #NOTE: sign?
    h_gas_wall = C_gas_wall * k_gas_wall_film * (X_gas_wall**n_gas_wall) * (T_tank_gas - T_wall_gas)

    Q_dot_vap_wall = h_gas_wall * A_gas_wall * (T_tank_gas - T_wall_gas)

    return Q_dot_vap_wall

#Q_dot_GS - convection
def solve_Q_dot_ullage_gas_liquid_surface():

    return 1

#Q_dot_VLC - condensation
def solve_Q_dot_condensation(T, m_dot_condensation):
    return m_dot_condensation*latent_heat_liq_vap(T)

#Q_dot_VLB - condensation
def solve_Q_dot_vaporization(T, m_dot_vaporization):
    return m_dot_vaporization*latent_heat_liq_vap(T)

#Q_dot_LW - convection
def solve_Q_dot_wall_liq():

    return 1

#Q_dot_LS - convection
def solve_Q_dot_liq_liq_surface_layer():

    return 1


#Q_dot_AWL
def solve_Q_env_liq_wall():

    return 1


#Q_dot_AWG
def solve_Q_env_vap_wall():

    return 1


def solve_m_dot_evaporated_liq():
    #NOTE: fluid restricted to one component and surface temperature assumed to equal saturation temperature
    h_ls = 
    
    E = 2.1e4 #correction coeff for N2O
    h_lsb = E*h_ls
    
    m_dot_evaporated_liq = h_lsb * (Area/latent_heat_liq_vap(T,MW))*(T_liq-T_sat)
    return m_dot_evaporated_liq


def solve_F_1(T, Z, A, B, alpha, n2o):

    #NOTE: n2o.MW --> g/mol | 
    F_1 = ((R_u*T)/n2o.MW) * (np.ln( (Z+2.414*B)/(Z-0.414*B) ) * (A/(5.657*B)) * (n2o.kappa/(n2o.Tc*alpha)) * (np.sqrt(alpha/T_ref) + n2o.kappa) )
    return F_1


def solve_F_2(T, Z, rho, A, B, alpha, n2o):

    F_2 = ((R_u*T)/n2o.MW) * ( (-1)*(Z/rho) * (A/((Z**2)+2*B*Z-(B**2))) * (1+n2o.kappa*np.sqrt(T_ref/alpha)) ) #NOTE: alpha from EOS?
    return F_2


def vap_phase_T_dot(T, rho, rho_dot, m, m_dot_prop, n2o, u_e, Q_dot_gas_wall):

    #convert rho to V_m
    V_m = (n2o.MW/1000) / rho

    #actually likely makes sense to keep T separate or define a new PR object????
    pr_eos_vap = PR(T=T, V=V_m, Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega) #TODO: check units
    P = pr_eos_vap.P #allegedly this is how pressure is solved (on init?)

    Z = (P*V_m)/(R_u*T)

    A = (pr_eos_vap.a*P)/( ((R_u/n2o.MW)**2)*T**2)
    B = (pr_eos_vap.b*P)/((R_u/n2o.MW)*T)

    alpha = (1+ n2o.kappa * (1 - np.sqrt(T/n2o.T_cr))) 

    F_1 = solve_F_1(T, Z, A, B, alpha, n2o)
    F_2 = solve_F_2(T, Z, A, B, alpha, rho, n2o)

    cv_ideal_gas = cp_ideal_gas(T) - (R_u/n2o.MW)

    T_dot =  (Q_dot_gas_wall - m_dot_prop*( (P/rho) + 0.5*(u_e**2) ) - m*F_2*rho_dot) / (m*(F_1+cv_ideal_gas))
    return T_dot


def liq_phase_T_dot_vap():

    T_dot = 1
    return T_dot


def liq_phase_T_dot_liq():

    T_dot = 1
    return T_dot



class model():

    #NOTE/TODO: THIS IS PROBABLY THE SAME as bens, change if not

    def __init__(self, oxidizer, TIMESTEP, T_amb, m_nos, Cd_1, A_inj_1, V_tank, Diam_tank, P_tank, P_cc, all_error, inj_model):
        self.oxidizer = oxidizer
        self.TIMESTEP = TIMESTEP #s
        self.T_amb = T_amb #K
        self.m_nos = m_nos #kg
        self.Cd_1 = Cd_1
        self.A_inj_1 = A_inj_1 #m^2
        self.V_tank = V_tank #m^3
        self.cross_sect_area = 0.25*np.pi*(Diam_tank**2)
        self.P_tank = P_tank #Pa
        self.all_error = all_error #%
        self.inj_model = inj_model

        self.n2o = Chemical('nitrous oxide')

        #setup - start by assuming the tank is in thermal equillibrium

        #solve molar volume
        V_m = (self.n2o.MW/1000) / (self.m_nos/self.V_tank)

        #input into PR EOS to solve temp :/

        #BUG: does this solve temperature?
        self.T_liq = PR(T=None, V=V_m, Tc=self.n2o.Tc, Pc=self.n2o.Pc, omega=self.n2o.omega)
        self.T_vap = self.T_liq #assuming tank starts in thermal equillibrium

        # Initialize nitrous oxide using the Chemical class to obtain critical properties
        self.pr_eos_liq = None #= PR(T=T, V=V_m, Tc=self.n2o.Tc, Pc=self.n2o.Pc, omega=self.n2o.omega)
        self.pr_eos_vap = None

        #TODO: solve quality of tank:

        self.m_liq = None
        self.m_vap = None

        self.rho_exit = 1



    def inst(self, P_cc):

        #TODO: PUT INJECTOR MODEL HERE!!!!
        #starting with spi model like the thesis used, once i prove everything runs i will swap in a better model here
        
        self.m_dot = spi_model(self.Cd_1, self.A_inj_1, self.P_tank, self.P_cc, self.rho_exit)



        if self.m_liq >= 0: #two phases in the tank!

            #assuming only liquid drains from the tank
            self.m_liq -= self.TIMESTEP * self.m_dot

            ###solve heat transfer terms



            ###solve T_dot 

            sol = solve_ivp( lambda T: liq_phase_T_dot_liq
















        else: #only gas in tank



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




