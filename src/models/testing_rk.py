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




def calculate_RHS( P:float, T:float, rho:float, m:float, P_dot_liq:float, T_dot:float, rho_dot_liq:float, m_dot:float, V_dot_liq:float ):

    F_1 = solve_F_1(T, rho)
    F_2 = solve_F_2(T, rho)

    cv_ig = (n2o.HeatCapacityGas.T_dependent_property(T)/MW) - (R_U/MW)

    u = CP.PropsSI('U', 'T', T, 'P', P, 'N2O') #J/kg

    V = m/rho #m^3

    RHS = (P*V_dot_liq + P_dot_liq*V) + m_dot*u + m*(F_1 + cv_ig)*T_dot + m*F_2*rho_dot_liq
    return RHS




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





def calculate_LHS(P_1_l, T_1_l, rho_1_l, m_dot_inj):

    h_l = CP.PropsSI('H', 'T', T_1_l, 'P', P_1_l, 'N2O')

    m_dot_evap = solve_m_dot_evaporated_liq(T_1_l, P_1_l, rho_1_l)
    Q_dot_evap = solve_Q_dot_evap(T_1_l, P_1_l, m_dot_evap)

    h_evap = solve_latent_heat_evap(T_1_l, P_1_l)

    Q_dot_ls = solve_Q_dot_ls(T_1_l, P_1_l, rho_1_l)


    #print(f"LHS: Q_dot_evap={Q_dot_evap} (J/s), Q_dot_ls={Q_dot_ls} (J/s), m_dot={m_dot} (kg/s), h_l={h_l} (J/kg), m_dot_evap={m_dot_evap}, h_evap={h_evap} (J/kg)\n")
    #print("checking LHS SIGN", (Q_dot_evap -Q_dot_ls), (m_dot*h_l), ( m_dot_evap*h_evap) )
    
    lhs = (Q_dot_evap -Q_dot_ls) + m_dot_inj*h_l + m_dot_evap*h_evap
    return lhs






def system_of_liq_odes(t, y, P_init, m_dot_inj, m_dot_evap, Q_dot_net, m, rho, rho_dot_liq):

    T, v = y  # Unpack state variables

    #print("before the crash: ", T, P_init, v)


    h_evap = (CP.PropsSI('H', 'T', T, 'D', rho, "N2O") - CP.PropsSI('H', 'T', T, 'D', rho, "N2O")) #J/kg

    V_dot_liq = solve_liq_phase_V_dot_liq( (v*m), rho, rho_dot_liq, m_dot_inj, m_dot_evap) 

    P_dot_liq = 0 #TODO: PASS THIS IN from previous P_dot_liq

    print("FIRST TIMESTEP T_DOT_LIQ: ", P_init, T, rho, Q_dot_net, m_dot_inj, m_dot_evap, h_evap, m, P_dot_liq, 0, rho_dot_liq)
    T_dot_liq = solve_liq_phase_T_dot_liq(P_init, T, rho, Q_dot_net, m_dot_inj, m_dot_evap, h_evap, m, P_dot_liq, 0, rho_dot_liq) #NOTE: just to use previous val V_dot_liq = 0
    
    print(V_dot_liq, T_dot_liq)
    return [T_dot_liq, V_dot_liq]





#### THIS WILL BECOME INIT
Cd_spi = 0.425
P_back = 1e5 #Pa
A_inj = 0.00003 #m^3 NOTE: GUESS

# Initial conditions NOTE: assuming tank starts at sat conditions
TIMESTEP = 1e-4
convergence_percent = 0.001 #NOTE: ADJUST AS REQUIRED,

P_init = 45e5 # Initial pressure in Pa
#T_init = 286.5
V_tank = 0.0354
m_nos = 20
v_init = V_tank/m_nos  # Initial volume in m^3

rho_sat_vap = CP.PropsSI('D', 'P', P_init, 'Q', 1, "N2O")
rho_sat_liq = CP.PropsSI('D', 'P', P_init, 'Q', 0, "N2O")
rho_sat_tank = 1/v_init

x_tank = ( (1/rho_sat_tank)-(1/rho_sat_liq) ) / ( (1/rho_sat_vap)-(1/rho_sat_liq) )

m_vap = x_tank*m_nos
m_liq = m_nos-m_vap

v_vap = 1/rho_sat_vap
V_vap = v_vap*m_vap
P_vap = P_init

v_liq = 1/rho_sat_liq
P_liq = P_init

rho_dot_liq = 0
P_dot_liq = 0
V_dot_liq = 0


Q_dot_evap = 0
m_dot_evap = 0
m_dot_inj = 1 
        
P_vap_prev = P_init
P_liq_prev = P_init

rho_vap_prev = 1/v_vap
rho_liq_prev = 1/v_liq

P_dot_vap = 0
rho_dot_vap = 0


#I THINK THIS IS JUST SETUP FOR TEST CASE BECAUSE THEY DIDNT MEASURE TEMPERATURE, MAY OR MAY NOT NEED IN ACTUAL SCRIPT INPUT
pr_eos = PR(P=P_liq, V= (MW*v_liq), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
T_liq = pr_eos.solve_T(P_liq, (MW*v_liq) )
T_vap = T_liq #starting at saturation

#### THIS WILL BECOME INST


### SOLVE THE LIQUID CV:

# Initial state setup
y0_liq = [T_liq, (v_liq*m_liq)]

# Solve the system

m_liq-=m_dot_inj*TIMESTEP

Q_dot_net_liq = Q_dot_evap #solve simple case to start
solution = solve_ivp(system_of_liq_odes, [0, TIMESTEP], y0_liq, args=(P_liq, m_dot_inj, m_dot_evap, Q_dot_net_liq, m_liq, rho_sat_liq, rho_dot_liq),
                     method='RK45', rtol=1e-8, atol=1e-8)

# Extract results
t = solution.t
y_liq = solution.y

T_liq = y_liq[0, -1]
V_liq = y_liq[1, -1]

pr_eos = PR(T=T_liq, V= (MW*V_liq/m_liq), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
P_liq = pr_eos.P

#print(f"t = {t[-1]:.5f}, T = {y[0, -1]:.4f}, V = {y[1, -1]}, P = {P_liq}, rho = {1/y[1, -1]}") # P = {y[1, i]:.4f},


# Print results (test code)
for i in range(len(t)):

    pr_eos = PR(T=y_liq[0, -1], V= (MW*y_liq[1, -1]/m_liq), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)

    print(f"t = {t[i]:.8f}, T = {y_liq[0, i]:.4f}, V = {y_liq[1, i]}, P = {P_liq}, rho = {1/y_liq[1, i]}") # P = {y[1, i]:.4f},





#### NOW RESOLVE LIQ TIME DERIVATIVES FOR VAPOR PHASE AND NEXT RUNGE KUTTA STEP

#TODO: Fill in, these will be stored as attributes

rho_liq = m_liq/V_liq

m_dot_evap = solve_m_dot_evaporated_liq(T_liq, P_liq, rho_liq)
Q_dot_evap = solve_Q_dot_evap(T_liq, P_liq, m_dot_evap)

#NOTE: P_vap? P_liq?
m_dot_inj = spi_model(Cd_spi, A_inj, P_vap, P_back, rho_liq) #TODO: FILL IN FOR TEST CASE THEY USED SPI AS WELL

rho_dot_liq = (rho_liq - rho_liq_prev)/TIMESTEP
V_dot_liq = solve_liq_phase_V_dot_liq(V_liq, rho_liq, rho_dot_liq, m_dot_inj, m_dot_evap)

h_evap = solve_latent_heat_evap(T_liq, P_liq)

Q_dot_net_liq = Q_dot_evap #simplified model just to get it to work at the start

#   huge error here vvvv


print("INTO NEXT TIMESTEP T_DOT_LIQ: ", P_liq, T_liq, rho_liq, Q_dot_net_liq, m_dot_inj, m_dot_evap, h_evap, m_liq, P_dot_liq, V_dot_liq, rho_dot_liq)
T_dot_liq = solve_liq_phase_T_dot_liq(P_liq, T_liq, rho_liq, Q_dot_net_liq, m_dot_inj, m_dot_evap, h_evap, m_liq, P_dot_liq, V_dot_liq, rho_dot_liq)
P_dot_liq = (P_liq - P_liq_prev)/TIMESTEP

print("P_dot_liq: ", P_dot_liq)

### VAPOR CV
V_dot_vap = (-V_dot_liq) #--> use this to solve rho_vap below!!!!!

#fwd euler step to solve rho_vap
print("checking sign convention: ", v_vap)
V_vap += V_dot_vap*TIMESTEP
print("checking sign convention after: ", v_vap)

Q_dot_net_vap = Q_dot_evap #simplified model just to get it to work at the start

solution_vap = solve_ivp(solve_liq_phase_T_dot_vap, [0, TIMESTEP], [T_vap], args=(P_vap, 1/(v_vap), Q_dot_net_vap, m_dot_inj, m_dot_evap, h_evap, m_vap, P_dot_vap, V_dot_vap, rho_dot_vap),
                     method='RK45', rtol=1e-8, atol=1e-8)

# Extract results
t = solution_vap.t
y_vap = solution_vap.y
T_vap = y_vap[0, -1]
V_vap = V_tank - V_liq

pr_eos = PR(T=T_vap, V= (MW*V_vap/m_vap), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
P_vap = pr_eos.P

###resolve derivatives
P_dot_vap = (P_vap - P_vap_prev)/TIMESTEP
#NOTE: already have V_dot_vap above 
rho_vap = m_vap/V_vap
rho_dot_vap = (rho_vap - rho_vap_prev)/TIMESTEP

T_dot_vap = solve_liq_phase_T_dot_vap(None, T_vap, P_vap, 1/(v_vap), Q_dot_net_vap, m_dot_inj, m_dot_evap, h_evap, m_vap, P_dot_vap, V_dot_vap, rho_dot_vap),





print("after", T_vap, P_vap, rho_vap, V_vap)

print("vap derivatives: ", T_dot_vap, V_dot_vap, rho_dot_vap, P_dot_vap)






#once we get here put inside criteria to switch to vapor phase