from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical
import matplotlib.pyplot as plt
import traceback

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

### solve latent heat of evaporation
"""
def solve_latent_heat_evap(T_gas, T_sat, T_liq, P_tank):

    #superheated gas
    if T_sat <= T_gas:


    #subcooled gas - metastable
    elif T_sat > T_gas: 

    return h_lv
"""




"""        ### solving U_dot_inj:
        n2o_ig = Chemical('N2O', T=T_liq) 
        preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
        T_sat = preos_l.Tsat(P_tank)
        #T_sat = CP.PropsSI('T', 'P', P_tank, 'Q', 1, 'N2O')
        h_liq = preos_l.H_dep_l/MW
        #latent_heat_evap_l = preos_l.Hvap(T_liq)/MW 

        preos_sat = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_sat, P=P_tank)
        h_sat_l = preos_sat.H_dep_l/MW #departure
        h_sat_g = preos_sat.H_dep_g/MW #departure

        preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
        #latent_heat_cond_g = (-1)*preos_g.Hvap(T_gas)/MW
        h_gas = preos_g.H_dep_g/MW #departure

        preos_sat = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_sat, P=P_tank)
        h_sat_l = preos_sat.H_dep_l/MW #+ n2o_ig.Cpg*(T_liq - T_REF)
        h_sat_g = preos_sat.H_dep_g/MW

        latent_heat_evap_l = preos_sat.Hvap(T_liq)/MW 
        #allegedly  ^ falls apart near the critical point
        delta_h_evap = ( (h_sat_g - h_sat_l) + (h_sat_l - h_liq) )"""


T_gas = 292.86796444582785
T_sat =  291.94857482143544
T_liq =  293.07156840787906
P_tank = 4938449.140297322


#currently in script!
n2o_ig = Chemical('N2O', T=T_liq) 
preos_l = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
T_sat = preos_l.Tsat(P_tank)
h_liq = preos_l.H_dep_l/MW

preos_sat = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_sat, P=P_tank)
h_sat_l = preos_sat.H_dep_l/MW #departure
h_sat_gas = preos_sat.H_dep_g/MW #departure

preos_g = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)
h_gas = preos_g.H_dep_g/MW #departure

lv_current = ( (h_sat_gas - h_sat_l) + (h_gas-h_liq) )

"""
###latent heat of vaporization!

#method 1:
#std but solved differently to see if there is anything unexpected
#h_lv = h_gas(T_sat) - h_liq(T_sat) since liquid surface at T_sat
preos_sat_1 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_sat, P=P_tank)
h_sat_l_dep = preos_sat_1.H_dep_l/MW
h_sat_g_dep = preos_sat_1.H_dep_g/MW

h_lv_1 = ( h_sat_g_dep - h_sat_l_dep )

#method 2:
#std
preos_sat_2 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_sat, P=P_tank)
h_lv_2 = preos_sat_2.Hvap(T_sat)/MW 
#allegedly  ^ falls apart near the critical point


#method 3:
#difference between enthalpies
n2o_ig_l_3 = Chemical('N2O', T=T_liq, P=P_tank) 
preos_l_3 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)

h_liq_3 = preos_l_3.H_dep_l/MW + ( n2o_ig_l_3.Cpg*(T_liq - T_REF) )

n2o_ig_g_3 = Chemical('N2O', T=T_gas, P=P_tank) 
preos_g_3 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_gas, P=P_tank)

h_gas_3 = preos_g_3.H_dep_g/MW + ( n2o_ig_g_3.Cpg*(T_gas - T_REF) )

h_lv_3 = ( h_gas_3 - h_liq_3 )

### solve latent heat of condensation
#def solve_latent_heat_cond():


print("current: ", lv_current, "1: ", h_lv_1, "2: ", h_lv_2, "3: ", h_lv_3)

"""

#### NIST inputs

R_U = 8.31446 #J/(mol K)
T_atm = 298.15 #K
P_atm = 101325 #Pa
T_crit = 309.52#K
P_crit = 7.2450
T_REF = T_crit
P_REF = P_crit
T_liq = 283
P_tank = 5e6





try:
    #current method: in sscript
    n2o_ig_l_1 = Chemical('N2O', T=T_liq) 
    preos_l_1 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    u_liq_1 = preos_l_1.U_dep_l/MW + n2o_ig_l_1.H - (R_U*(T_liq)/MW)
except Exception as e:
    u_liq_1 = 0

try:
    #method 2: including T_REF
    n2o_ig_l_2 = Chemical('N2O', T=T_liq) 
    preos_l_2 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    u_liq_2 = preos_l_2.U_dep_l/MW + ( n2o_ig_l_2.H - (R_U*(T_liq-T_REF)/MW) )
except Exception as e:
    u_liq_2 = 0

try:
    #method 3: try difference of ideal gas enthalpy as well?
    n2o_ig_l_3 = Chemical('N2O', T=T_liq) 
    n2o_ig_l_ref_3 = Chemical('N2O', T=T_REF) 
    preos_l_3 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    u_liq_3 = preos_l_3.U_dep_l/MW + ( (n2o_ig_l_3.H- n2o_ig_l_3.H )- (R_U*(T_liq-T_REF)/MW) )
except Exception as e:
    u_liq_3 = 0

#method 4: this makes more sense than 3? 
try:
    n2o_ig_l_ref_4 = Chemical('N2O', T=T_REF) 
    preos_l_4 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_REF, P=P_REF)
    u_4_ref = preos_l_4.U_dep_l/MW + n2o_ig_l_ref_4.H - (R_U*(T_REF)/MW)

    n2o_ig_l_4 = Chemical('N2O', T=T_liq) 
    preos_l_4 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_REF)
    u_4 = preos_l_4.U_dep_l/MW + n2o_ig_l_ref_4.H - (R_U*(T_liq)/MW)

    u_liq_4 = u_4 - u_4_ref
except Exception as e:
    u_liq_4 = 0

try:
    #method 5: mixing up signs?
    n2o_ig_l_ref_5 = Chemical('N2O', T=T_REF) 
    preos_l_5 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_REF, P=P_REF)
    u_5_ref = preos_l_5.U_dep_l/MW + n2o_ig_l_ref_5.H - (R_U*(T_REF)/MW)

    n2o_ig_l_5 = Chemical('N2O', T=T_liq) 
    preos_l_5 = PR(Tc=TC, Pc=PC, omega=OMEGA, T=T_liq, P=P_tank)
    u_5 = preos_l_5.U_dep_l/MW + n2o_ig_l_ref_5.H - (R_U*(T_liq)/MW)

    u_liq_5 = u_5_ref - u_5
except Exception as e:
    u_liq_5 = 0

print("testing internal energy: ", u_liq_1, u_liq_2, u_liq_3, u_liq_4, u_liq_5 )
u_NIST_WEBBOOK = 179.84 *1000 #J/kg
print("expecting: ", u_NIST_WEBBOOK )

#NOTE: FROM NIST WEBBOOK ENTHALPY CONVENTION H=0 AT the normal boiling point T = 184.68