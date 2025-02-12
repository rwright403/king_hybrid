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
u_gas_arr = []
u_liq_arr = []
for T in T_arr:
    n2o_ig = Chemical('N2O', T=T) 

    try:
        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T, P=P_tank)
        u_gas = preos_g.U_dep_g/MW + (n2o_ig.H - (R_U*T/MW) )
        u_gas_arr.append(u_gas)
    except Exception as e:
        u_gas_arr.append(0)

    try:
        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T, P=P_tank)
        u_liq = preos_g.U_dep_l/MW + (n2o_ig.H - (R_U*T/MW) )
        u_liq_arr.append(u_liq)
    except Exception as e:
        u_liq_arr.append(0)




plt.scatter(T_arr, u_gas_arr, label = 'gas' )
plt.scatter(T_arr, u_liq_arr, label = 'liq' )
plt.xlabel('Temp (K)')
plt.ylabel('Internal Energy (J)')
plt.title('Internal Energy vs. Temp')
plt.legend()
plt.grid(True)
plt.show()