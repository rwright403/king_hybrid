from scipy.optimize import minimize
from thermo import Chemical
from thermo import IdealGas
from thermo import PR
import numpy as np
from rocketprops import RocketProps
import CoolProp as CP

# Global Constants:
R_U = 8.31446 #J/(mol K) 

T_REF = 298.15 #K
P_REF = 101325 #Pa
rp_n2o = RocketProps("N2O") # for viscosity only
n2o = Chemical('nitrous oxide')
n2o_ideal = IdealGas('N2O')

MW = (n2o.MW/1000) #n2o.MW in g/mol --> convert to kg/mol
KAPPA = 0.37464 + 1.5422*n2o.omega - 0.26992*n2o.omega**2

TANK_DIAM = 0.0254*5.5 #m
CS_AREA = 0.25*np.pi*(TANK_DIAM**2) #m^2
g = 9.81 #m/s^2



# Heat Transfer Functions and LHS
def latent_heat_vap(T): #TODO: convert to coolprop
    return n2o.H_gas(T) - n2o.H_liquids(T)

def solve_Q_dot_evap(T, m_dot_evap):
    return m_dot_evap*latent_heat_vap(T)

def solve_h_ls(T_1_l, P_1_l, rho_1_l, T_1_s):
    K_H = ?
    C = ?
    n = ?

    k_l = CP.PropsSI('k', 'T', T_1_l, 'P', P_1_l, 'N2O')
    Cp_l = CP.PropsSI('CPMASS', 'T', T_1_l, 'P', P_1_l, 'N2O')
    visc_l = rp_n2o.viscosity(T_1_l, P_1_l) #iirc cp doesnt have this property so i used rocketprops

    dV_dT_P = n2o.get_partial_derivative('V', 'T', T=T_1_l, P=P_1_l)
    beta = dV_dT_P*rho_1_l

    Gr =  ((TANK_DIAM**3)*(rho_1_l**2)*g*beta*np.abs(T_1_l - T_1_s) ) / (visc_l**2)
    Pr = (Cp_l*visc_l)/ k_l

    X_ls = Gr*Pr
    h_ls = K_H * C * (k_l/TANK_DIAM) * X_ls**n
    return h_ls

#Q_dot_LS - convection between bulk liquid surface and liquid surface layer
def solve_Q_dot_ls(T_1_l, P_1_l, rho_1_l):

    T_1_s = CP.PropsSI('T', 'P', P_1_l, 'Q', 0, 'N2O')
    h_ls = solve_h_ls(T_1_l, P_1_l, rho_1_l, T_1_s)

    Q_dot_ls = h_ls*CS_AREA*(T_1_l - T_1_s) #liquid temp minus surface temp --> assuming surface temp is sat temp?
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


def calculate_LHS(T_1_l, P_1_l, rho_1_l, m_dot_inj, h_l):

    m_dot_evap = solve_m_dot_evaporated_liq(T_1_l, P_1_l, rho_1_l)
    Q_dot_evap = solve_Q_dot_evap(T_1_l, m_dot_evap)

    h_evap = latent_heat_vap(T_1_l)

    Q_dot_ls = solve_Q_dot_ls(T_1_l, P_1_l, rho_1_l)

    lhs = (-Q_dot_evap -Q_dot_ls) - m_dot_inj*h_l - m_dot_evap*h_evap
    return lhs

# Thermodynamic Functions and RHS

def solve_F_1(T, Z, A, B, alpha):

    F_1 = ((R_U*T)/MW) * (np.log( (Z+2.414*B)/(Z-0.414*B) ) * (A/(5.657*B)) * (n2o.KAPPA/(n2o.Tc*alpha)) * (np.sqrt(alpha/T_REF) + KAPPA) )
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
    T_dot = (T-T_1)/TIMESTEP

    Z = (P*V_m)/(R_U*T)

    #does calling PR EOS for this still work with constraints?
    alpha = (1+ KAPPA * (1 - np.sqrt(T/n2o.Tc))) 

    a = (0.4572 * ((R_U/MW)**2) * (n2o.Tc**2)/n2o.Pc) * alpha
    b = 0.07780 * (R_U/MW) * n2o.Tc/n2o.Pc

    A = (a*P)/( ((R_U/n2o.MW)**2)*T**2)
    B = (b*P)/((R_U/n2o.MW)*T)

    F_1 = solve_F_1(T, Z, A, B, alpha, n2o)
    F_2 = solve_F_2(T, Z, A, B, alpha, rho, n2o)

    cv_ig = n2o_ideal.cp(T) - (R_U/MW)

    #NEED TO SOLVE u for liquid here!!!!! coolprop?
    u = CP.PropsSI('U', 'T', T, 'P', P, 'N2O')

    RHS = (P*V_dot + P_dot*(rho*m)) + m_dot*u + m*(F_1 + cv_ig)*T_dot + m*F_2*rho_dot
    return RHS


# Optimization Problem

# Initial Setup of Previous Values:
P_1 = 
T_1 = 
rho_1 =
m =
m_dot =
TIMESTEP =

time_constants = (P_1, T_1, rho_1, m, m_dot, TIMESTEP)

# Step 1: Define the objective function
def objective(x):
    P, T, rho = x
    LHS = calculate_LHS(P, T, rho)  # Your function to compute LHS
    RHS = calculate_RHS(P, T, rho)  # Your function to compute RHS
    return (LHS - RHS)**2  # Minimizing the square of the absolute difference


# Step 2: Define phase constraints (all inequality constraints)
def phase_constraint_p(x):
    P = x[0]
    return n2o.Pc - P  # Ensures P <= P_crit

def phase_constraint_t(x):
    T = x[1]
    return n2o.Tc - T  # Ensures T <= T_crit

def phase_constraint_rho(x):
    rho = x[2]
    return rho - n2o.rhoc  # Ensures rho >= rho_crit


# Step 3: Define initial guesses and constraints list

T_guess = .999*T_1
rho_guess = .999*rho_1

#solving pressure w pr eos to only feed optimizer physically possible vals
pr_eos = PR(T=T_guess, V_m=(MW/ rho_guess), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
print(pr_eos.P) #i think this works like this
P_guess = pr_eos.P

initial_guess = [P_guess, T_guess, rho_guess] 

constraints = [
    {'type': 'ineq', 'fun': phase_constraint_p},      # P <= P_crit
    {'type': 'ineq', 'fun': phase_constraint_t},      # T <= T_crit
    {'type': 'ineq', 'fun': phase_constraint_rho}     # rho >= rho_crit
]


# Step 4: Run the optimizer using SLSQP
result = minimize(objective, initial_guess, args=time_constants, constraints=constraints, method='SLSQP')

# Analyze the results
if result.success:
    P_opt, T_opt, rho_opt = result.x
    print(f"Optimal values: P={P_opt}, T={T_opt}, rho={rho_opt}")
else:
    print("Optimization failed:", result.message)