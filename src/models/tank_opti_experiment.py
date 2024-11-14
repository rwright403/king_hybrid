from scipy.optimize import minimize
from thermo import Chemical
from thermo import PR
import numpy as np
import rocketprops
import CoolProp.CoolProp as CP

# Global Constants:
R_U = 8.31446 #J/(mol K) 

T_REF = 298.15 #K
P_REF = 101325 #Pa
n2o = Chemical('nitrous oxide')

MW = (n2o.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o.omega - 0.26992*n2o.omega**2

TANK_DIAM = 0.0254*5.5 #m
CS_AREA = 0.25*np.pi*(TANK_DIAM**2) #m^2
g = 9.81 #m/s^2


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
def latent_heat_vap(T): #TODO: convert to coolprop
    return (CP.PropsSI('H', 'T', T, 'Q', 1, "N2O") - CP.PropsSI('H', 'T', T, 'Q', 0, "N2O")) #J/kg

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


def calculate_LHS(T_1_l, P_1_l, rho_1_l, m_dot):

    h_l = CP.PropsSI('H', 'T', T_1_l, 'P', P_1_l, 'N2O')

    m_dot_evap = solve_m_dot_evaporated_liq(T_1_l, P_1_l, rho_1_l)
    Q_dot_evap = solve_Q_dot_evap(T_1_l, m_dot_evap)

    h_evap = latent_heat_vap(T_1_l)

    Q_dot_ls = solve_Q_dot_ls(T_1_l, P_1_l, rho_1_l)


    #print(f"LHS: Q_dot_evap={Q_dot_evap} (J/s), Q_dot_ls={Q_dot_ls} (J/s), m_dot={m_dot} (kg/s), h_l={h_l} (J/kg), m_dot_evap={m_dot_evap}, h_evap={h_evap} (J/kg)\n")
    #print("checking LHS SIGN", (Q_dot_evap -Q_dot_ls), (m_dot*h_l), ( m_dot_evap*h_evap) )
    
    lhs = (Q_dot_evap -Q_dot_ls) + m_dot*h_l + m_dot_evap*h_evap
    return lhs

# Thermodynamic Functions and RHS

def solve_F_1(T, Z, A, B, alpha):

    F_1 = ((R_U*T)/MW) * (np.log( (Z+2.414*B)/(Z-0.414*B) ) * (A/(5.657*B)) * (KAPPA/(n2o.Tc*alpha)) * (np.sqrt(alpha/T_REF) + KAPPA) )
    return F_1


def solve_F_2(T, Z, rho, A, B, alpha):

    F_2 = ((R_U*T)/MW) * ( (-1)*(Z/rho) * (A/((Z**2)+2*B*Z-(B**2))) * (1+KAPPA*np.sqrt(T_REF/alpha)) )
    return F_2


def calculate_RHS(P, T, rho, P_1, T_1, rho_1, m, m_dot, TIMESTEP):

    #TODO: what is m_dot and u???? i think u would be internal energy at prev step

    V_m = MW/ rho

    V = m / rho
    V_1 = m / rho_1


    #estimating derivatives from previous values
    V_dot = (V - V_1)/TIMESTEP
    rho_dot = (rho - rho_1)/TIMESTEP
    P_dot = (P-P_1)/TIMESTEP
    #print("checking pressure deriv: ", P, P_1, (P-P_1), P_dot)
    T_dot = (T-T_1)/TIMESTEP

    Z = (P*V_m)/(R_U*T)

    #does calling PR EOS for this still work with constraints?
    alpha = (1+ KAPPA * (1 - np.sqrt(T/n2o.Tc))) 

    a = (0.4572 * ((R_U/MW)**2) * (n2o.Tc**2)/n2o.Pc) * alpha
    b = 0.07780 * (R_U/MW) * n2o.Tc/n2o.Pc

    A = (a*P)/( ((R_U/MW)**2)*T**2)
    B = (b*P)/((R_U/MW)*T)
    #print(b, P, (R_U/MW), T)

    F_1 = solve_F_1(T, Z, A, B, alpha)
    F_2 = solve_F_2(T, Z, A, B, alpha, rho)

    #n2o_ideal = IdealGas('N2O', T=T_1_l)
    #print("CP", n2o.HeatCapacityGas.T_dependent_property(T_1_l)/MW)
    cv_ig = (n2o.HeatCapacityGas.T_dependent_property(T)/MW) - (R_U/MW)
    #BUG: i think the error is here


    #NEED TO SOLVE u for liquid here!!!!! coolprop?
    u = CP.PropsSI('U', 'T', T, 'P', P, 'N2O') #J/kg

    #print(f"RHS: P={P}, T={T}, rho={rho}")

    print(f"RHS: P={P}, V_dot={V_dot}, P_dot={P_dot}, rho={rho}, m={m},m_dot={m_dot}, u={u}, F_1={F_1}, cv_ig={cv_ig}, T_dot={T_dot}, F_2={F_2}, rho_dot={rho_dot}")
    print("RHS", P*V_dot, P_dot*(m/rho), m_dot*u, m*(F_1 + cv_ig)*T_dot, m*F_2*rho_dot)

    #print(m_dot*u, m_dot, u)

    RHS = (P*V_dot + P_dot*(m/rho)) + m_dot*u + m*(F_1 + cv_ig)*T_dot + m*F_2*rho_dot
    return RHS


# Optimization Problem

# Initial Setup of Previous Values:
T_1_l = 300
rho_1_l = 700

#need to solve this wiht PR EOS
pr_eos = PR(T=T_1_l, V=(MW/rho_1_l), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
P_1_l = pr_eos.P

# print(f"INPUTS: P={P_1_l} (Pa), T={T_1_l} (K), rho={rho_1_l}, (kg/m^3)")

m = 1
m_dot = -1 #kg/s
TIMESTEP = 0.05 #s

time_constants = (P_1_l, T_1_l, rho_1_l, m, m_dot, TIMESTEP)

# Step 1: Define the objective function
def objective(x, *args):
    T, rho = x
    pr_eos = PR(T=T_guess, V=(MW/rho_guess), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
    P = pr_eos.P
    #print("in obj function", T,rho,P)

    # Unpack constants
    P_1_l, T_1_l, rho_1_l, m, m_dot, TIMESTEP = args

    LHS = calculate_LHS(T_1_l, P_1_l, rho_1_l, m_dot) # Your function to compute LHS
    RHS = calculate_RHS(P, T, rho, P_1_l, T_1_l, rho_1_l, m, m_dot, TIMESTEP)  # Your function to compute RHS
    print(f"\nOBJECTIVE: {(LHS - RHS)**2}, LHS={LHS}, RHS={RHS}\n")
    return (LHS - RHS) # (LHS - RHS)**2  # Minimizing the square of the absolute difference


# Step 2: Define phase constraints
"""
def phase_constraint_T(x):
    T = x[0]
    return n2o.Tc - T  # Ensures T <= T_crit

def phase_constraint_rho(x):
    rho = x[1]
    return rho - n2o.rhoc  # Ensures rho >= rho_crit
"""
# Define your bounds for T and rho
bounds = [(182.33, n2o.Tc),  # T bounded between 0 and T_crit
          (n2o.rhoc, None)]  # rho bounded between rho_crit and infinity (or some large upper limit)

# Step 3: Define initial guesses and constraints list

T_guess = .95*T_1_l
rho_guess = .95*rho_1_l

#solving pressure w pr eos to only feed optimizer physically possible vals
pr_eos = PR(T=T_guess, V=(MW/rho_guess), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)

initial_guess = [T_guess, rho_guess] 

"""
constraints = [
    {'type': 'ineq', 'fun': phase_constraint_T},      # T <= T_crit
    {'type': 'ineq', 'fun': phase_constraint_rho}     # rho >= rho_crit
]
"""

# Step 4: Run the optimizer using SLSQP
#result = minimize(objective, initial_guess, args=time_constants, constraints=constraints, method='SLSQP')
result = minimize(objective, initial_guess, args=time_constants, bounds=bounds, method='Powell')

# Analyze the results
if result.success:
    T_opt, rho_opt = result.x

    #solve P under PR EOS
    pr_eos = PR(T=T_opt, V=(MW/rho_opt), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
    P_opt = pr_eos.P

    print(f"Optimal values: P={P_opt}, T={T_opt}, rho={rho_opt}")
else:
    print("Optimization failed:", result.message)