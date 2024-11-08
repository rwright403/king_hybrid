#I love coolprop! ~ units: http://www.coolprop.org/v4/apidoc/CoolProp.html

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from rocketprops.rocket_prop import get_prop
import numpy as np
from thermo import Chemical
from thermo.eos import PR


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
    
def compressibility(A, B):
    # Define the coefficients of the polynomial in terms of A and B
    coeffs = [1, -(1 - B)**2, (A - 3 * B**2 - 2 * B), -(A * B - B**2 - 3 * B)]
    
    # Calculate the roots
    roots = np.roots(coeffs)
    
    # Only return the largest real root
    real_roots = [root.real for root in roots if root.imag == 0]
    
    return max(real_roots) if real_roots else None

def cp_sat_vap(T):
    MW = 44.013  # g/mol
    RU = 8.314  # J/(mol·K)
    RM = RU / (MW / 1000)  # J/(kg·K)
    T_crit = 309.56  # K
    P_crit = 7.238e6  # Pa
    omega = 0.160  # Acentric factor (unitless)

    T_ref = T / T_crit  # Unitless

    kappa = 0.37464 + 1.54226 * omega - 0.26992 * (omega ** 2)  # Unitless

    alpha = (1 + kappa * (1 - np.sqrt(T_ref))) ** 2  # Unitless
    a_crit = (0.45724 * (RU ** 2) * (T_crit ** 2)) / P_crit  # J·m^3 / mol^2

    a = a_crit * alpha  # J·m^3 / mol^2
    b = 0.07780 * RU * T_crit / P_crit  # m^3 / mol

    P_sat = P_sat_vap(T) #Pa
    A = a * P_sat / ((RM**2)*(T**2))  # Unitless
    B = b * P_sat / (RM * T)  # Unitless

    Z = compressibility(A, B)  # Unitless

    pr_eos_vap = PR(T=T, P=P_sat, Tc=T_crit, Pc=P_crit, omega=omega)
    rho = (MW/1000) / pr_eos_vap.V_g  # kg/m^3

    F_2 = (RM * T) * ((-Z / rho) * (A / ((Z**2) + 2 * B * Z - (B**2))) * (1 + kappa * np.sqrt(T_ref / alpha)))  # J/kg

    print(F_2, rho, P_sat, A, B, kappa, alpha, RM, T)
    
    #setup
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
        T/=1000
        val = 1000*( A + B*T + C*(T**2) + D*(T**3) + E/(T**2) ) / MW #J/(kg K)
        return F_2*val #J/(kg K)
    

pObj = get_prop('N2O')

Temps = np.linspace(188,300,50)

cp_arr = []
polynomial_arr = []
rocketprops_arr = []
for T in Temps:
    polynomial_arr.append(cp_sat_vap(T))
    cp_arr.append(CP.PropsSI('CPMASS', 'Q', 0, 'T', T, 'N2O'))
    rocketprops_arr.append(4186.8*pObj.CpAtTdegR(1.8*T)) #convert to rankine for rocketprops input

plt.plot(Temps, cp_arr, label = 'COOLPROP')
plt.plot(Temps, polynomial_arr, label = 'POLYNOMIAL')
plt.plot(Temps, rocketprops_arr, label = 'ROCKETPROPS')
plt.xlabel('temp (K)')
plt.ylabel('CPMASS (J/(kg K)')
plt.title('Thermodynamic Property Comparison')
plt.grid(True)
plt.legend()
plt.show()