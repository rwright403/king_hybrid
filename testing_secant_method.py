import numpy as np
import matplotlib.pyplot as plt

# Define the function for TWOPHASEerror
def TWOPHASEerror(eta_crit, omega):
    function_diff = (eta_crit**2) + ((omega**2)-2*omega)*((1-eta_crit)**2) + 2*(omega**2)*np.log(eta_crit) + 2*(omega**2)*(1-eta_crit)  
    return function_diff

# Define the secant method
def secant(func, x1):
    x_eps = x1 * 0.005  # Set the tolerance to be 0.5% of init guess
    x2 = x1 -x1 * 0.01  # Set a second point 1% away from the original guess
    F1 = func(x1)  # Evaluate function at x1
    F2 = func(x2)  # Evaluate function at x2
    kk = 1  # Set up counter
    kk_max = 1000


    while np.abs(x2 - x1) >= (x_eps) and kk < kk_max:  # While error is too large and counter is less than max
        x3 = x2 - (F2 * (x2 - x1) / (F2 - F1)) 
        x1 = x2  # Move everything forward
        x2 = x3
        F1 = F2
        F2 = func(x2)
        kk = kk + 1
    x = x2
    return x

# Example parameters
omega_sat = 3.3746557566463062 # Replace with your actual value
eta_crit_sat = 0.9328827548961321  # Initial guess for critical pressure ratio
all_err = 1e-4  # Error tolerance

# Step 1: Loop to find the critical pressure ratio
while np.abs(TWOPHASEerror(eta_crit_sat, omega_sat) ) > all_err:
        eta_crit_sat = secant((lambda T: TWOPHASEerror(T, omega_sat)), eta_crit_sat)

# Output the result
print(f'Critical pressure ratio found: eta_crit_sat = {eta_crit_sat:.5f}')

# Plotting the TWOPHASEerror function
eta_values = np.linspace(0.01, 1.0, 500)  # Range of eta_crit values from 0.01 to 1.0
error_values = [TWOPHASEerror(eta, omega_sat) for eta in eta_values]

plt.figure(figsize=(10, 6))
plt.plot(eta_values, error_values, label='TWOPHASEerror', color='blue')
plt.axhline(0, color='red', linestyle='--', label='y=0 (Root)')
plt.axvline(eta_crit_sat, color='green', linestyle='--', label='Estimated Root')
plt.title('TWOPHASEerror Function')
plt.xlabel('eta_crit')
plt.ylabel('TWOPHASEerror Value')
plt.ylim([-10, 10])  # Set limits for better visualization
plt.xlim([0, 1])
plt.grid()
plt.legend()
plt.show()