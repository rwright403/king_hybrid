from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical
import matplotlib.pyplot as plt
import traceback

T_REF = 298.15 #K
P_REF = 1e5 #Pa
R_U = 8.314 #J/mol/K

n2o = Chemical('nitrous oxide', T=T_REF)

MW = (n2o.MW/1000) #n2o.MW in g/mol --> converted to kg/mol

#NOTE: T_min = 298 K, T_max = 1400 K
#https://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Mask=1

def solve_cp_g_nist(T):
    A = 27.67988
    B = 51.14898
    C = -30.64454
    D = 6.847911
    E = -0.157906
    F = 71.24934
    G = 238.6164
    H = 82.04824
    t = T/1000

    cp = A +B*t + C*t**2 + D*t**3 + E/t**2
    return cp/MW #convert to J/kg K

#NOTE: This returns enthalpy relative to T_REF = 298.15 and P_REF = 1 atm.
def solve_h_g_nist(T):
    A = 27.67988
    B = 51.14898
    C = -30.64454
    D = 6.847911
    E = -0.157906
    F = 71.24934
    #G = 238.6164
    H = 82.04824
    t = T/1000

    h = A*t + (B/2)*t**2 + (C/3)*t**3 + (D/4)*t**4 - E/t + F - H
    return (h/MW)*1000 #convert to J/kg from kJ/mol

T_1 = 287.944 #K
P_1 = 4.5e6 #Pa

#NOTE: testing gas n2o

preos_c = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_1, P=P_1)
print("06-jan-2025 thermo property test")
print(preos_c.phase)




n2o_ig_ref = Chemical('N2O', T=T_REF) 
preos_ref = PR(Tc=n2o_ig_ref.Tc, Pc=n2o_ig_ref.Pc, omega=n2o_ig_ref.omega, T=T_REF, P=P_REF)

cp_ref_ig = n2o_ig_ref.Cpg

n2o_ig_1 = Chemical('N2O', T=T_1) 
preos_1 = PR(Tc=n2o_ig_1.Tc, Pc=n2o_ig_1.Pc, omega=n2o_ig_1.omega, T=T_1, P=P_1)

cp_1_ig = n2o_ig_1.Cpg

h_test_1 = cp_1_ig*(T_1) - cp_ref_ig*T_REF #+ preos_1.H_dep_g/MW
#h_nist = solve_h_g_nist(T_1)
#h_test_nist = solve_cp_g_nist(T_1)*T_1 - solve_cp_g_nist(T_REF)*T_REF

h_test_2 = n2o_ig_1.H
#print(n2o_ig_1.U)

#print(h_test_1, h_nist, h_test_nist, h_test_2)
#print( n2o_ig_1.Cpg, solve_cp_g_nist(T_1) )
#print("negative enthalpy test: ", cp_1_ig*T_1+preos_1.H_dep_g/MW) 

preos_ref =  PR(Tc=n2o_ig_1.Tc, Pc=n2o_ig_1.Pc, omega=n2o_ig_1.omega, T=T_REF, P=P_REF)

h_real_test_1 = preos_1.H_dep_l/MW + n2o_ig_1.H 
h_real_test_2 = preos_1.H_dep_l/MW + n2o_ig_1.H  + (R_U*T_REF)/MW



h_real_cp = CP.PropsSI('H', 'P', P_1, 'T', T_1, 'N2O')
cool_state = CP.PhaseSI('P', P_1, 'T', T_1, 'N2O')
print("coolprop state: ", cool_state)

#print(h_real_test_1, h_real_test_2, h_real_cp )


u_ig = (n2o_ig_1.H - (R_U*T_1/MW) )
u_real_test = preos_1.U_dep_l/MW + u_ig

u_cool = CP.PropsSI('U', 'P', P_1, 'T', T_1, 'N2O') - CP.PropsSI('U', 'T', T_REF, 'P', P_REF, 'N2O')

print(u_real_test, u_cool)

preos_l = PR(Tc=n2o_ig_1.Tc, Pc=n2o_ig_1.Pc, omega=n2o_ig_1.omega, T=T_1, P=P_1)
partial_du_dv_const_T_liq = T_1*preos_l.dP_dT_l - P_1 #NOTE: derivative is partial at constant volume
n2o_ig_l = Chemical('N2O', T=T_1) 
    
u_liq = preos_l.U_dep_l/MW + (n2o_ig_l.H - (R_U*T_1/MW) )#NOTE: same warning as ^ but also for this line

print("model? ", u_liq)

