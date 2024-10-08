import numpy as np
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt

#testing hem model fix

#inputs: (grab from real vals in script)
P_tank = 5147690.721566
P_cc = 296662.0
h_tank_exit = 244162.19860526677
h_inj_exit = 199138.2128330915
Cd_1 = 0.6
A_inj_1 = 3.090605599220321e-5

m_dot_ox = None

s_inj = CP.PropsSI('S', 'H', h_tank_exit, 'P', P_tank, 'N2O')
h_inj_exit = CP.PropsSI('H', 'S', s_inj, 'P', P_cc, 'N2O')
            
#rho_exit = CP.PropsSI('D', 'S', s_inj, 'P', P_cc, 'N2O')
#m_dot_ox  = Cd_1 * A_inj_1 * rho_exit * np.sqrt( 2 * (h_tank_exit -  h_inj_exit) )


downstream_pres_arr = np.linspace(1e5, P_tank, 100)
m_dot_HEM_arr = []

for pres in downstream_pres_arr:
    rho_exit = CP.PropsSI('D', 'S', s_inj, 'P', pres, 'N2O')
    m_dot_HEM = Cd_1 * A_inj_1 * rho_exit * np.sqrt( 2 * (h_tank_exit -  h_inj_exit) )
    m_dot_HEM_arr.append(m_dot_HEM)

m_dot_HEM_crit = np.max(m_dot_HEM_arr)
P_crit = downstream_pres_arr[np.argmax(m_dot_HEM_arr)]

if P_cc < P_crit:
    print("choked flow")
    m_dot_ox = m_dot_HEM_crit

else:
    print("unchoked")
    rho_exit = CP.PropsSI('D', 'S', s_inj, 'P', P_cc, 'N2O')
    m_dot_ox = Cd_1 * A_inj_1 * rho_exit * np.sqrt( 2 * (h_tank_exit -  h_inj_exit) )

#print(m_dot_ox, rho_exit, "choked rho: ", CP.PropsSI('D', 'S', s_inj, 'P', P_crit, 'N2O'), "unchoked rho: ", CP.PropsSI('D', 'S', s_inj, 'P', P_cc, 'N2O'))
#print("h --> f(P_cc): ", h_inj_exit, "h --> f(P_crit): ", CP.PropsSI('H', 'S', s_inj, 'P', P_crit, 'N2O') )

"""
delta_p_arr = []
for p in downstream_pres_arr:
    delta_p_arr.append(P_tank-p)


plt.plot(delta_p_arr,m_dot_HEM_arr)
plt.xlabel('delta P (MPa)')
plt.ylabel('Mass Flow Rate (kg/s)')
plt.title('delta P (MPa) vs Mass Flow Rate (kg/s)')
plt.grid(True)
plt.show()
"""