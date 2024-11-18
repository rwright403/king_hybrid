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






def liq_phase_V_dot_liq(V, rho, rho_dot, m_dot_inj, m_dot_evap):

    V_dot = (-m_dot_inj - m_dot_evap - V*rho_dot) / rho #NOTE: NOT SPECIFIC VOLUME!!!!
    print("V_dot", V_dot)
    return V_dot


    
def liq_phase_T_dot_liq(P, T, rho, Q_dot_net, m_dot_inj, m_dot_evap, h_evap, m, P_dot, V_dot, rho_dot):

    F_1 = solve_F_1(T, rho)
    F_2 = solve_F_2(T, rho)

    u = CP.PropsSI('U', 'T', T, 'P', rho, 'N2O') #J/kg
    h = CP.PropsSI('U', 'T', T, 'P', P, 'N2O') #J/kg

    cv_ideal_gas = (n2o.HeatCapacityGas.T_dependent_property(T)/MW) - (R_U/MW)

    #NOTE: FIRST TIMESTEP STARTS IN EQUILIBRIUM
    V_dot = 0
    P_dot = 0
    rho_dot = 0

    V = m/rho

    T_dot = (Q_dot_net - m_dot_inj*h - m_dot_evap*h_evap -(P*V_dot + P_dot*V) - m_dot_inj*u - m*F_2*rho_dot ) / (m*(F_1+cv_ideal_gas))

    print("T_dot", T_dot)
    return T_dot





def calculate_RHS( P:float, T:float, rho:float, m:float, P_dot:float, T_dot:float, rho_dot:float, m_dot:float, V_dot:float ):

    F_1 = solve_F_1(T, rho)
    F_2 = solve_F_2(T, rho)

    cv_ig = (n2o.HeatCapacityGas.T_dependent_property(T)/MW) - (R_U/MW)

    u = CP.PropsSI('U', 'T', T, 'P', P, 'N2O') #J/kg

    V = m/rho #m^3

    RHS = (P*V_dot + P_dot*V) + m_dot*u + m*(F_1 + cv_ig)*T_dot + m*F_2*rho_dot
    return RHS





def latent_heat_vap(T,P): 
    return (CP.PropsSI('H', 'T', T, 'P', P, "N2O") - CP.PropsSI('H', 'T', T, 'P', P, "N2O")) #J/kg
#BUG: assuming equilibrium here, that seems wrong... #NOTE: YES BAD !!!

def solve_Q_dot_evap(T, P, m_dot_evap):
    #print("m_dot_evap", m_dot_evap)
    return m_dot_evap*latent_heat_vap(T, P) #J/s

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
    m_dot_evaporated_liq = h_lsb * (CS_AREA/latent_heat_vap(T_1_l, P_1_l))*(T_1_l-T_1_s)
    return m_dot_evaporated_liq





def calculate_LHS(P_1_l, T_1_l, rho_1_l, m_dot_inj):

    h_l = CP.PropsSI('H', 'T', T_1_l, 'P', P_1_l, 'N2O')

    m_dot_evap = solve_m_dot_evaporated_liq(T_1_l, P_1_l, rho_1_l)
    Q_dot_evap = solve_Q_dot_evap(T_1_l, P_1_l, m_dot_evap)

    h_evap = latent_heat_vap(T_1_l, P_1_l)

    Q_dot_ls = solve_Q_dot_ls(T_1_l, P_1_l, rho_1_l)


    #print(f"LHS: Q_dot_evap={Q_dot_evap} (J/s), Q_dot_ls={Q_dot_ls} (J/s), m_dot={m_dot} (kg/s), h_l={h_l} (J/kg), m_dot_evap={m_dot_evap}, h_evap={h_evap} (J/kg)\n")
    #print("checking LHS SIGN", (Q_dot_evap -Q_dot_ls), (m_dot*h_l), ( m_dot_evap*h_evap) )
    
    lhs = (Q_dot_evap -Q_dot_ls) + m_dot_inj*h_l + m_dot_evap*h_evap
    return lhs






def system_of_odes(t, y, P_init, m_dot_inj, m_dot_evap, Q_dot_net, m, rho, rho_dot):
    """
    System of ODEs for T, P, V.
    
    y: State variables [T, P, V]
    """
    T, v = y  # Unpack state variables

    #print("before the crash: ", T, P_init, v)


    h_evap = (CP.PropsSI('H', 'T', T, 'D', rho, "N2O") - CP.PropsSI('H', 'T', T, 'D', rho, "N2O")) #J/kg

    V_dot = liq_phase_V_dot_liq( (v*m), rho, rho_dot, m_dot_inj, m_dot_evap) 

    P_dot = 0 #TODO: PASS THIS IN from previous P_dot
    
    T_dot = liq_phase_T_dot_liq(P_init, T, rho, Q_dot_net, m_dot_inj, m_dot_evap, h_evap, m, P_dot, 0, rho_dot)
    
    print(V_dot, T_dot)
    return [T_dot, V_dot]







# Initial conditions NOTE: assuming tank starts at sat conditions
TIMESTEP = 1e-4
convergence_percent = 0.001 #NOTE: ADJUST AS REQUIRED,

P_init = 45e5 # Initial pressure in Pa
#T_init = 286.5
Tank_V = 0.0354
m_nos = 20
v_init = Tank_V/m_nos  # Initial volume in m^3



rho_sat_vap = CP.PropsSI('D', 'P', P_init, 'Q', 1, "N2O")
rho_sat_liq = CP.PropsSI('D', 'P', P_init, 'Q', 0, "N2O")
rho_sat_tank = 1/v_init

x_tank = ( (1/rho_sat_tank)-(1/rho_sat_liq) ) / ( (1/rho_sat_vap)-(1/rho_sat_liq) )

m_vap = x_tank*m_nos
m_liq = m_nos-m_vap

v_vap = 1/rho_sat_vap
P_vap = P_init
    

v_liq = 1/rho_sat_liq
P_liq = P_init

rho_dot = 0
P_dot = 0
V_dot = 0


Q_dot_evap = 0
m_dot_evap = 0
m_dot_inj = 1 #NOTE: i think this sign is true, if it isnt the script is broken and it makes it appear true
h_evap = (CP.PropsSI('H', 'P', P_init, 'Q', 1, "N2O") - CP.PropsSI('H', 'P', P_init, 'Q', 0, "N2O")) #J/kg
        
P_vap_prev = P_init
P_liq_prev = P_init



pr_eos = PR(P=P_init, V= (MW*v_liq), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
T_init = pr_eos.solve_T(P_init, (MW*v_liq) )

print(T_init, P_init, v_liq)

# Initial state
y0 = [T_init, (v_liq*m_liq)]

# Solve the system
print(m_liq)
m_liq-=m_dot_inj*TIMESTEP
print(m_liq)

solution = solve_ivp(system_of_odes, [0, TIMESTEP], y0, args=(P_init, m_dot_inj, m_dot_evap, 0, m_liq, rho_sat_liq, rho_dot),
                     method='RK45', rtol=1e-8, atol=1e-8)

# Extract results
t = solution.t
y = solution.y
pr_eos = PR(T=y[0, -1], V= (MW*y[1, -1]/m_liq), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
P_liq = pr_eos.P
T_liq = y[0, -1]
V_liq = y[1, -1]

print(f"t = {t[-1]:.5f}, T = {y[0, -1]:.4f}, V = {y[1, -1]}, P = {P_liq}, rho = {1/y[1, -1]}") # P = {y[1, i]:.4f},


# Print results
for i in range(len(t)):

    pr_eos = PR(T=y[0, -1], V= (MW*y[1, -1]/m_liq), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)

    print(f"t = {t[i]:.8f}, T = {y[0, i]:.4f}, V = {y[1, i]}, P = {P_liq}, rho = {1/y[1, i]}") # P = {y[1, i]:.4f},



#### NOW RESOLVE TIME DERIVATIVES FOR NEXT TIMESTEP


