import numpy as np
import matplotlib.pyplot as plt
from thermo import Chemical
from thermo.eos import PR

# Global Constants:
R_U = 8.31446 #J/(mol K)
T_REF = 298.15 #K

n2o_g = Chemical('nitrous oxide', T=T_REF)
PC = n2o_g.Pc
TC = n2o_g.Tc
OMEGA = n2o_g.omega

MW = (n2o_g.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o_g.omega - 0.26992*n2o_g.omega**2

P_tank = 2.5e6 #Pa

T_lower_bound = 270 #K
T_upper_bound = 310 #K

T_arr = np.linspace(T_lower_bound, T_upper_bound, 50)
cv_gas_arr = []
cv_liq_arr = []
for T in T_arr:
    n2o_ig = Chemical('N2O', T=T) 

    try:
        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T, P=P_tank)
        Cv_gas = preos_g.Cv_dep_g/MW + n2o_g.Cvg
        cv_gas_arr.append(Cv_gas)
    except Exception as e:
        cv_gas_arr.append(0)

    try:
        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T, P=P_tank)
        Cv_liq = preos_l.Cv_dep_l/MW + n2o_g.Cvg
        cv_liq_arr.append(Cv_liq)
    except Exception as e:
        cv_liq_arr.append(0)




plt.scatter(T_arr, cv_gas_arr, label = 'gas' )
plt.scatter(T_arr, cv_liq_arr, label = 'liq' )
plt.xlabel('Temp (K)')
plt.ylabel('Specific Heat at const vol (J/ (kg K))')
plt.title('Specific Heat vs. Temp')
plt.legend()
plt.grid(True)
plt.show()