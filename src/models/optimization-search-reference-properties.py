import numpy as np
from scipy.optimize import minimize
from thermo import Chemical  # Assuming you're using the thermo library for fluid properties
from thermo.eos import PR

# Fluid
fluid = 'N2O'

# Global Constants:
R_U = 8.31446 #J/(mol K)
T_REF = 298.15 #K

n2o_global = Chemical('nitrous oxide', T=T_REF)
PC = n2o_global.Pc
TC = n2o_global.Tc
OMEGA = n2o_global.omega

MW = (n2o_global.MW/1000) #n2o.MW in global/mol --> converted to kglobal/mol
KAPPA = 0.37464 + 1.5422*n2o_global.omega - 0.26992*n2o_global.omega**2
b = 0.07780*(R_U*TC/PC)
g = 9.81 #m/s^2


# Define the thermodynamic property function (e.g., enthalpy or internal energy)
# This should return the property based on T and P
def calculate_enthalpy(T, P):
    preos = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T, P=P)

    return preos.H_dep_g

# Objective function to minimize (difference from 0 for reference enthalpy)
def objective_function(x):
    T, P = x  # T: Temperature (K), P: Pressure (atm)
    
    # Calculate the enthalpy at the current T and P
    enthalpy = calculate_enthalpy(T, P)
    
    # We want enthalpy to be as close to 0 as possible (reference point)
    return abs(enthalpy)

# Optimization bounds for temperature and pressure (adjust as needed)
bounds = [(1, 600), (1e3, 3e6)]  # Temperature in K, Pressure in Pa

# Initial guess for optimization (midpoint of bounds)
initial_guess = [250, 1e6]  # T = 250 K, P = 1 atm (this is just an initial guess)

# Run the optimization to minimize the objective function
result = minimize(objective_function, initial_guess, bounds=bounds, method='L-BFGS-B')

# Print the result
T_opt, P_opt = result.x
optimal_enthalpy = calculate_enthalpy(T_opt, P_opt)

print(f"Optimal Reference Temperature: {T_opt:.2f} K")
print(f"Optimal Reference Pressure: {P_opt:.2f} Pa")
print(f"Enthalpy at this point: {optimal_enthalpy:.2f} J/kg")

