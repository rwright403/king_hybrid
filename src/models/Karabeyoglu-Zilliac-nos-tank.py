from thermo import Chemical
from thermo.eos import PR

# Create a Fluid object for nitrous oxide
nitrous_oxide = Chemical('nitrous oxide')

from scipy.integrate import solve_ivp
import numpy as np

#global constants:
R_u = 8.31446 #J/K/mol


def spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit):
    m_dot_spi = Cd_hem_spi_dyer * A_inj_ox * np.sqrt( 2 * rho_tank_exit * (P_1 - P_2)  )
    #print(rho_tank_exit)
    return m_dot_spi

#property relations

def cp_ideal_gas(T):
    A = 21.62
    B = 78.81
    C = -57.78
    D = 18.3
    E = 0

    T_max = 310 #K
    T_min = 150 #K

    if not (T_min <= T <= T_max):
        raise ValueError(f"Input T={T} is out of bounds! Cp (ideal gas) relation only valid for [{T_min},{T_max}].")
    else:
        cp_ideal_gas = A + B*T + C*(T**2) + D*(T**3) + E*(T**4)
        return cp_ideal_gas



def solve_F_1(T, Z, n2o):
    A = 1
    B = 1
    #NOTE: n2o.MW --> g/mol | 

    F_1 = ((R_u*T)/n2o.MW) * (np.ln( (Z+2.414*B)/(Z-0.414*B) ) * (A/(5.657*B)) * (n2o.kappa/(n2o.Tc*alpha)) * (np.sqrt(alpha/T_ref) + n2o.kappa) )
    return F_1


def solve_F_2(T, Z, rho, n2o):
    A = 1
    B = 1
    F_2 = ((R_u*T)/n2o.MW) * ( (-1)*(Z/rho) * (A/((Z**2)+2*B*Z-(B**2))) * (1+n2o.kappa*np.sqrt(T_ref/alpha)) ) #NOTE: alpha from EOS?
    return F_2

def gas_phase_T_dot(T, rho, rho_dot, m, m_dot_prop, n2o, u_e, Q_dot_gas_wall):

    #convert rho to V_m
    V_m = n2o.MW / rho

    #actually likely makes sense to keep T separate or define a new PR object????
    pr_eos_vap = PR(T=T, V=V_m, Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega) #TODO: check units
    P = pr_eos_vap.P #allegedly this is how pressure is solved (on init?)

    Z = (P*V_m)/(R_u*T)

    F_1 = solve_F_1(T, Z, n2o)
    F_2 = solve_F_2(T, Z, rho, n2o)

    cv_ideal_gas = cp_ideal_gas() - (R_u/n2o.MW)

    T_dot =  (Q_dot_gas_wall - m_dot_prop*( (P/rho) + 0.5*(u_e**2) ) - m*F_2*rho_dot) / (m*(F_1+cv_ideal_gas))
    return T_dot

def liq_phase_T_dot_gas():

    T_dot = 1
    return T_dot


def liq_phase_T_dot_liq():

    T_dot = 1
    return T_dot



class model():

    #NOTE/TODO: THIS IS PROBABLY THE SAME as bens, change if not

    def __init__(self, oxidizer, TIMESTEP, m_liq, Cd_1, A_inj_1, V_tank, P_tank, P_cc, all_error, inj_model):
        self.oxidizer = oxidizer
        self.TIMESTEP = TIMESTEP
        self.Cd_1 = Cd_1
        self.A_inj_1 = A_inj_1
        self.V_tank = V_tank
        self.P_tank = P_tank
        self.P_cc = P_cc
        self.all_error = all_error
        self.inj_model = inj_model

        #setup

        # Initialize nitrous oxide using the Chemical class to obtain critical properties
        self.n2o = Chemical('nitrous oxide')
        self.pr_eos_liq = None #= PR(T=T, V=V_m, Tc=self.n2o.Tc, Pc=self.n2o.Pc, omega=self.n2o.omega)
        self.pr_eos_vap = None

        #TODO: solve quality of tank:
        self.m_liq = None
        self.m_vap = None

        self.m_tank = self.m_liq+self.m_vap

        self.rho_exit = 1



    def inst(self, P_cc):

        #TODO: PUT INJECTOR MODEL HERE!!!!
        #starting with spi model like the thesis used, once i prove everything runs i will swap in a better model here
        
        self.m_dot = spi_model(self.Cd_1, self.A_inj_1, self.P_tank, self.P_cc, self.rho_exit)



        if self.m_liq >= 0: #two phases in the tank!
            print("todo: implement")




        else: #only gas in tank



            #update conservation of mass

            self.m_tank -= self.TIMESTEP * self.m_dot

            rho_gas_prev = self.rho_gas
            self.rho_gas = self.m_tank/self.V_tank

            rho_dot = (self.rho_gas - rho_gas_prev) / self.TIMESTEP
            

            """ --> not sure how this works but Runge-Kutta to solve T here!
            # Solve the system using Runge-Kutta
            sol = solve_ivp(lambda T: , TIMESTEP, y0, method='RK45', t_eval=t_eval)

            # Extract results
            T = sol.y[0]  # Temperature over time

            """

            sol = solve_ivp(lambda T: , self.TIMESTEP, y0, method='RK45', t_eval=t_eval)

            #now that we solved T and have rho, can update tank pressure with Peng Robinson EOS

            #peng robinson eos takes molar volume, convert rho to molar volume

            # Initialize the PR EOS with the given temperature and molar volume
            self.pr_eos_vap = PR(T=T, V=V_m, Tc=self.n2o.Tc, Pc=self.n2o.Pc, omega=self.n2o.omega)
            P = self.pr_eos_vap.P
            #resolve all thermodynamic properties











        

