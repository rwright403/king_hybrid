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






def liq_phase_v_dot_liq(V, rho, rho_dot, m_dot_inj, m_dot_evap, m):
    v_dot = (-m_dot_inj - m_dot_evap - V*m*rho_dot) / (m*rho) #NOTE: NOT SPECIFIC VOLUME!!!!
    print("v_dot",v_dot)
    return v_dot

"""
def liq_phase_P_dot_liq(P, P_prev, TIMESTEP): #NOTE: I DONT LIKE THIS!   - Bad
    P_dot_liq = (P - P_prev) / TIMESTEP

    return P_dot_liq
"""
    
def liq_phase_T_dot_liq(P, T, rho, Q_dot_net, m_dot_inj, m_dot_evap, h_evap, m):

    F_1 = solve_F_1(T, rho)
    F_2 = solve_F_2(T, rho)

    u = CP.PropsSI('U', 'T', T, 'P', rho, 'N2O') #J/kg
    h = CP.PropsSI('U', 'T', T, 'P', P, 'N2O') #J/kg

    cv_ideal_gas = (n2o.HeatCapacityGas.T_dependent_property(T)/MW) - (R_U/MW)

    #NOTE: FIRST TIMESTEP STARTS IN EQUILIBRIUM
    V_dot = 0
    P_dot = 0
    rho_dot = 0


    T_dot = (Q_dot_net - m_dot_inj*h - m_dot_evap*h_evap -(P*V_dot + P_dot/rho) - m_dot_inj*u - m*F_2*rho_dot ) / (m*(F_1+cv_ideal_gas))

    print("T_dot",T_dot)
    return T_dot












def system_of_odes(t, y, P_init, m_dot_inj, m_dot_evap, Q_dot_net, m, rho):
    """
    System of ODEs for T, P, V.
    
    y: State variables [T, P, V]
    """
    T, V = y  # Unpack state variables

    print("before the crash: ", T, P_init, V)


    h_evap = (CP.PropsSI('H', 'T', T, 'Q', 1, "N2O") - CP.PropsSI('H', 'T', T, 'Q', 0, "N2O")) #J/kg

    # Get the rate of change of temperature, pressure, and volume
    V_dot = liq_phase_v_dot_liq(V, rho, 0, m_dot_inj, m_dot_evap, m) #NOTE: RHO_dot zero for testing
    #rho_dot = (-1/(V**2))*V_dot

    #print("V_dot:" , V_dot)
    T_dot = liq_phase_T_dot_liq(P_init, T, rho, Q_dot_net, m_dot_inj, m_dot_evap, h_evap, m) 
    
    #pr_eos = PR(T=T, V=(MW*V), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
    #P = pr_eos.P
    #P_dot = liq_phase_P_dot_liq(P, P_prev, TIMESTEP)
    

    return [T_dot, V_dot]







# Initial conditions
P_init = 45e5 # Initial pressure in Pa
#T_init = 286.5
Tank_V = 0.0354
m_nos = 20
V_init = Tank_V/m_nos  # Initial volume in m^3



rho_sat_vap = CP.PropsSI('D', 'P', P_init, 'Q', 1, "N2O")
rho_sat_liq = CP.PropsSI('D', 'P', P_init, 'Q', 0, "N2O")
rho_sat_tank = m_nos/V_init

x_tank = ( (1/rho_sat_tank)-(1/rho_sat_liq) ) / ( (1/rho_sat_vap)-(1/rho_sat_liq) )

m_vap = x_tank*m_nos
m_liq = m_nos-m_vap

v_vap = 1/rho_sat_vap
P_vap = P_init
    

v_liq = 1/rho_sat_liq
P_liq = P_init


Q_dot_evap = 0
m_dot_evap = 0
m_dot_inj = 1 #NOTE: i think this sign is true, if it isnt the script is broken and it makes it appear true
h_evap = (CP.PropsSI('H', 'P', P_init, 'Q', 1, "N2O") - CP.PropsSI('H', 'P', P_init, 'Q', 0, "N2O")) #J/kg
        
P_vap_prev = P_init
P_liq_prev = P_init

# Time span for the integration
#t_end = 10
#t_eval = np.linspace(t0, t_end, 100)
TIMESTEP = 1e-4

pr_eos = PR(P=P_init, V= (MW*v_liq), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
T_init = pr_eos.solve_T(P_init,(MW*v_liq) )

# Initial state
y0 = [T_init, v_liq]

# Solve the system
#solution = solve_ivp(system_of_odes, [t0, t_end], y0, args=(m_dot_inj, m_dot_evap, Q_dot_net, m, rho, F_1, cv_ideal_gas, F_2, u, P_prev, TIMESTEP), method='RK45', t_eval=t_eval)

                                                            #P_init, m_dot_inj, m_dot_evap, Q_dot_net, m, rho)
print(m_liq)
m_liq-=m_dot_inj*TIMESTEP
print(m_liq)

solution = solve_ivp(system_of_odes, [0, TIMESTEP], y0, args=(P_init, m_dot_inj, m_dot_evap, 0, m_liq, rho_sat_liq),
                     method='RK45', rtol=1e-8, atol=1e-8)

# Extract results
t = solution.t
y = solution.y

# Print results
for i in range(len(t)):

    
    pr_eos = PR(T=y[0, i], V= (MW*y[1, i]), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)

    #print( v_liq, y[1, i])


    print(f"t = {t[i]:.8f}, T = {y[0, i]:.4f}, V = {y[1, i]}, P = {pr_eos.P}, rho = {1/y[1, i]}") # P = {y[1, i]:.4f},

