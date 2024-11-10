from scipy.optimize import minimize
from thermo import Chemical
from thermo import IdealGas
from thermo import PR
import numpy as np

# Global Constants:
R_U = 8.31446 #J/(mol K) 

T_REF = 298.15 #K
P_REF = 101325 #Pa

n2o = Chemical('nitrous oxide')
n2o_ideal = IdealGas('N2O')




# Heat Transfer Functions and LHS

def solve_Q_dot_evap(T, m_dot_evap):
    #latent heat from ideal gas relation under model?
    return m_dot_evap*latent_heat_liq_vap(T)


#Q_dot_LS - convection
def solve_Q_dot_ls():

    return 1

def solve_m_dot_evaporated_liq():
    #NOTE: fluid restricted to one component and surface temperature assumed to equal saturation temperature
    h_ls = 
    
    E = 2.1e4 #correction coeff for N2O from thesis
    h_lsb = E*h_ls
    
    #latent heat from ideal gas relation under model?
    m_dot_evaporated_liq = h_lsb * (Area/latent_heat_liq_vap(T,MW))*(T_liq-T_sat)
    return m_dot_evaporated_liq


def calculate_LHS(m_dot_inj,h_l):

    Q_dot_evap
    Q_dot_ls = solve_Q_dot_ls()

    m_dot_evap

    lhs = (-Q_dot_evap -Q_dot_ls) - m_dot_inj*h_l - m_dot_evap*h_evap
    return lhs

# Thermodynamic Functions and RHS

def solve_F_1(T, Z, A, B, alpha):

    #NOTE: n2o.MW --> g/mol | 
    F_1 = ((R_U*T)/(n2o.MW/1000)) * (np.log( (Z+2.414*B)/(Z-0.414*B) ) * (A/(5.657*B)) * (n2o.kappa/(n2o.Tc*alpha)) * (np.sqrt(alpha/T_REF) + n2o.kappa) )
    return F_1


def solve_F_2(T, Z, rho, A, B, alpha):

    F_2 = ((R_U*T)/n2o.MW) * ( (-1)*(Z/rho) * (A/((Z**2)+2*B*Z-(B**2))) * (1+n2o.kappa*np.sqrt(T_REF/alpha)) ) #NOTE: alpha from EOS?
    return F_2


def calculate_RHS(P, T, rho):

    #TODO: need to bring in V_1, P_1, T_1, m,
    #TODO: what is m_dot and u????

    #convert rho to V_m
    V_m = (n2o.MW/1000) / rho
    V_m_1 = (n2o.MW/1000) / rho_1

    #estimating derivatives from previous values
    V_dot = (V_m - V_1)/TIMESTEP
    rho_dot = (rho - rho_1)/TIMESTEP
    P_dot = (P-P_1)/TIMESTEP
    T_dot = (T-T_1)/TIMESTEP

    Z = (P*V_m)/(R_U*T)

    #does calling PR EOS for this still work with constraints?

    kappa = 0.37464 + 1.5422*n2o.omega - 0.26992*n2o.omega**2
    alpha = (1+ kappa * (1 - np.sqrt(T/n2o.Tc))) 

    a = (0.4572 * ((R_U/(n2o.MW/1000))**2) * (n2o.Tc**2)/n2o.Pc) * alpha
    b = 0.07780 * (R_U/(n2o.MW/1000)) * n2o.Tc/n2o.Pc

    A = (a*P)/( ((R_U/n2o.MW)**2)*T**2)
    B = (b*P)/((R_U/n2o.MW)*T)

    F_1 = solve_F_1(T, Z, A, B, alpha, n2o)
    F_2 = solve_F_2(T, Z, A, B, alpha, rho, n2o)

    cv_ig = n2o_ideal.cp(T)  - (R_U/(n2o.MW/1000))

    RHS = (P*V_dot + P_dot*(rho*m)) + m_dot*u + m*(F_1 + cv_ig)*T_dot + m*F_2*rho_dot
    return RHS


# Optimization Problem

# Step 1: Define the objective function
def objective(x):
    P, T, rho = x
    LHS = calculate_LHS(P, T, rho)  # Your function to compute LHS
    RHS = calculate_RHS(P, T, rho)  # Your function to compute RHS
    return (LHS - RHS)**2  # Minimizing the square of the absolute difference

# Step 2: Define the Peng-Robinson EOS constraint (equality constraint)
def peng_robinson_constraint(x):
    P, T, rho = x
    

    #TODO: IMPLEMENT

    return 1

# Step 3: Define phase constraints (inequality constraints)
def phase_constraint_p(x):
    P = x[0]
    return n2o.Pc - P  # Ensures P <= P_crit

def phase_constraint_t(x):
    T = x[1]
    return n2o.Tc - T  # Ensures T <= T_crit

def phase_constraint_rho(x):
    rho = x[2]
    rho_c = n2o.MW/(1000*n2o.Vc)
    return rho - rho_c  # Ensures rho >= rho_crit


# Initial Setup:







# Step 4: Define initial guesses and constraints list
initial_guess = [P_guess, T_guess, rho_guess]
constraints = [
    {'type': 'eq', 'fun': peng_robinson_constraint},  # Peng-Robinson EOS constraint
    {'type': 'ineq', 'fun': phase_constraint_p},      # P <= P_crit
    {'type': 'ineq', 'fun': phase_constraint_t},      # T <= T_crit
    {'type': 'ineq', 'fun': phase_constraint_rho}     # rho >= rho_crit
]

# Step 5: Run the optimizer using SLSQP
result = minimize(objective, initial_guess, constraints=constraints, method='SLSQP')

# Analyze the results
if result.success:
    P_opt, T_opt, rho_opt = result.x
    print(f"Optimal values: P={P_opt}, T={T_opt}, rho={rho_opt}")
else:
    print("Optimization failed:", result.message)



