import numpy as np
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt

#testing hem model fix

#inputs: (grab from real vals in script)
T_tank = 283
P_tank = 2.55e6+CP.PropsSI('P', 'Q', 0, 'T', T_tank, 'N2O')
P_cc_arr = np.linspace(296662.0, P_tank, 100)
h_tank_exit = 244162.19860526677
h_inj_exit = 199138.2128330915
Cd_1 = 0.66
A_inj_1 = 0.25*np.pi*((1.5e-3)**2)

m_dot_final_dyer_arr = []
m_dot_final_hem_arr = []
m_dot_final_spi_arr = []
delta_p_arr = []

for P_cc in P_cc_arr:

    ### start model

    # SPI MODEL
    rho_1_spi = CP.PropsSI('D', 'H', h_tank_exit, 'P', P_tank, 'N2O') 
    m_dot_spi = Cd_1 * A_inj_1 * np.sqrt( 2 * rho_1_spi * (P_tank - P_cc)  )

    # HEM MODEL
    m_dot_hem = None
    downstream_pres_arr = np.linspace(1e5, P_tank, 100)
    m_dot_HEM_arr = []

    for pres in downstream_pres_arr:
        s_inj = CP.PropsSI('S', 'H', h_tank_exit, 'P', P_tank, 'N2O')
        h_inj_exit = CP.PropsSI('H', 'S', s_inj, 'P', pres, 'N2O')
        rho_exit = CP.PropsSI('D', 'S', s_inj, 'P', pres, 'N2O')
        m_dot_HEM = Cd_1 * A_inj_1 * rho_exit * np.sqrt( 2 * (h_tank_exit -  h_inj_exit) )
        m_dot_HEM_arr.append(m_dot_HEM)
        #print(m_dot_HEM)


    m_dot_HEM_crit = np.max(m_dot_HEM_arr)
    P_crit = downstream_pres_arr[np.argmax(m_dot_HEM_arr)]

    

    if P_cc < P_crit:
        print("choked flow")
        m_dot_hem = m_dot_HEM_crit

        P_sat = CP.PropsSI('P', 'Q', 0, 'T', T_tank, 'N2O')
        rho_1_spi = CP.PropsSI('D', 'H', h_tank_exit, 'P', P_tank, 'N2O') 
        m_dot_spi = Cd_1 * A_inj_1 * np.sqrt( 2 * rho_1_spi * (P_tank - P_crit)  )
        dyer_k = np.sqrt( (P_tank - P_cc) / (P_sat - P_cc) ) 

        m_dot_dyer = ((dyer_k/(1+dyer_k)) * m_dot_spi) + ((1/(1+dyer_k)) * m_dot_hem)

    else:
        print("unchoked")
        s_inj = CP.PropsSI('S', 'H', h_tank_exit, 'P', P_tank, 'N2O')
        h_inj_exit = CP.PropsSI('H', 'S', s_inj, 'P', P_cc, 'N2O')
        rho_exit = CP.PropsSI('D', 'S', s_inj, 'P', P_cc, 'N2O')
        m_dot_hem = Cd_1 * A_inj_1 * rho_exit * np.sqrt( 2 * (h_tank_exit -  h_inj_exit) )

    P_sat = CP.PropsSI('P', 'Q', 0, 'T', T_tank, 'N2O')
    dyer_k = np.sqrt( (P_tank - P_cc) / (P_sat - P_cc) ) 

    print(dyer_k)


    m_dot_dyer = ((dyer_k/(1+dyer_k)) * m_dot_spi) + ((1/(1+dyer_k)) * m_dot_hem)
    
    
    
    m_dot_final_dyer_arr.append(m_dot_dyer)
    m_dot_final_hem_arr.append(m_dot_hem)
    m_dot_final_spi_arr.append(m_dot_spi)
    delta_p_arr.append(P_tank-P_cc)


plt.plot(delta_p_arr,m_dot_final_dyer_arr, color = 'rebeccapurple', label = 'dyer')
plt.plot(delta_p_arr,m_dot_final_hem_arr, color = 'r', label = 'hem')
plt.plot(delta_p_arr,m_dot_final_spi_arr, color = 'b', label = 'spi')

plt.xlabel('delta P (MPa)')
plt.ylabel('Mass Flow Rate (kg/s)')
plt.title('delta P (MPa) vs Mass Flow Rate (kg/s)')
plt.grid(True)
plt.show()



