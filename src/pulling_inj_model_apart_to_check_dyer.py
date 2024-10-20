#thank you for the incredible thesis: https://emersonvn.com/project/two_phase_injector/

#I love coolprop! ~ units: http://www.coolprop.org/v4/apidoc/CoolProp.html

import CoolProp.CoolProp as CP
from rocketprops.rocket_prop import get_prop #NOTE: just using because CP doesn't have nitrous viscosity
from rocketprops.rocket_prop import Propellant
import matplotlib.pyplot as plt
import numpy as np

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

def TWOPHASEerror(eta_crit, omega):
    function_diff = (eta_crit**2) + ((omega**2)-2*omega)*((1-eta_crit)**2) + 2*(omega**2)*np.log(eta_crit) + 2*(omega**2)*(1-eta_crit)  
    return function_diff

def LOWSUBCOOLEDerror(eta_crit, eta_sat, omega_sat):
    function_diff = (((omega_sat+(1/omega_sat)-2)/(2*eta_sat))*(eta_crit**2)) - (2*(omega_sat-1)*eta_crit) + (omega_sat*eta_sat*np.log(eta_crit/eta_sat)) + ((3/2)*omega_sat*eta_sat) - 1
    return function_diff

def thermo_span_wagner(rho, T, param):
    # Constants for N2O
    R = 8.3144598 / 44.0128 * 1000  # Gas constant (J/kg*K)
    T_c = 309.52  # Critical Temperature (K)
    rho_c = 452.0115  # Critical Density (kg/m^3)

    n0 = np.array([0.88045, -2.4235, 0.38237, 0.068917, 0.00020367, 0.13122, 0.46032,
          -0.0036985, -0.23263, -0.00042859, -0.042810, -0.023038])
    n1 = n0[0:5]
    n2 = n0[5:12]
    a1 = 10.7927224829
    a2 = -8.2418318753
    c0 = 3.5
    v0 = np.array([2.1769, 1.6145, 0.48393])
    u0 = np.array([879, 2372, 5447])
    t0 = np.array([0.25, 1.125, 1.5, 0.25, 0.875, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5])
    d0 = np.array([1, 1, 1, 3, 7, 1, 2, 5, 1, 1, 4, 2])
    P0 = np.array([1, 1, 1, 2, 2, 2, 3])
    t1 = t0[0:5]
    t2 = t0[5:12]
    d1 = d0[0:5]
    d2 = d0[5:12]

    # Calculate non-dimensional variables
    tau = T_c / T
    delta = rho / rho_c

    # Calculate explicit Helmholtz energy and derivatives
    ao = a1 + a2 * tau + np.log(delta) + (c0 - 1) * np.log(tau) + np.sum(v0 * np.log(1 - np.exp(-u0 * tau / T_c)))
    ar = np.sum(n1 * tau**t1 * delta**d1) + np.sum(n2 * tau**t2 * delta**d2 * np.exp(-delta**P0))
    ao_tau = a2 + (c0 - 1) / tau + np.sum(v0 * u0 / T_c * np.exp(-u0 * tau / T_c) / (1 - np.exp(-u0 * tau / T_c)))
    ao_tautau = -(c0 - 1) / tau**2 + np.sum(-v0 * u0**2 / T_c**2 * np.exp(-u0 * tau / T_c) / (1 - np.exp(-u0 * tau / T_c))**2)
    ar_tau = np.sum(n1 * t1 * tau**(t1 - 1) * delta**d1) + np.sum(n2 * t2 * tau**(t2 - 1) * delta**d2 * np.exp(-delta**P0))
    ar_tautau = np.sum(n1 * t1 * (t1 - 1) * tau**(t1 - 2) * delta**d1) + np.sum(n2 * t2 * (t2 - 2) * tau**(t2 - 2) * delta**d2 * np.exp(-delta**P0))
    ar_delta = np.sum(n1 * d1 * delta**(d1 - 1) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    ar_deltadelta = np.sum(n1 * d1 * (d1 - 1) * delta**(d1 - 2) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 2) * ((d2 - P0 * delta**P0) * (d2 - 1 - P0 * delta**P0) - P0**2 * delta**P0) * np.exp(-delta**P0))
    ar_deltatau = np.sum(n1 * d1 * t1 * delta**(d1 - 1) * tau**(t1 - 1)) + np.sum(n2 * t2 * tau**(t2 - 1) * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))

    out = 0.0
    if param == 'p':  # Pressure (Pa)
        out = rho * R * T * (1 + delta * ar_delta)
    elif param == 'u':  # Specific internal energy (J/kg)
        out = R * T * tau * (ao_tau + ar_tau)
    elif param == 's':  # Specific entropy (J/kg*K)
        out = R * (tau * (ao_tau + ar_tau) - ao - ar)
    elif param == 'h':  # Specific enthalpy (J/kg)
        out = R * T * (1 + tau * (ao_tau + ar_tau) + delta * ar_delta)
    elif param == 'cv':  # Specific heat constant volume (J/kg*K)
        out = R * -tau**2 * (ao_tautau + ar_tautau)
    elif param == 'cp':  # Specific heat constant pressure (J/kg*K)
        out = R * (-tau**2 * (ao_tautau + ar_tautau) + (1 + delta * ar_delta - delta * tau * ar_deltatau)**2 / (1 + 2 * delta * ar_delta + delta**2 * ar_deltadelta))
    elif param == 'a':  # Speed of sound (m/s)
        out = np.sqrt(R * T * (1 + 2 * delta * ar_delta + delta**2 * ar_deltadelta - (1 + delta * ar_delta - delta * tau * ar_deltatau)**2 / (tau**2 * (ao_tautau + ar_tautau))))
    else:
        raise ValueError('Invalid input')

    return out

def dyer_model( Cd_hem_spi_dyer, A_inj_ox, P_1, P_sat, P_2, h_1 ):

    #print("dyer check inputs (P1,P_sat,P_2)",P_1, P_sat, P_2, )

    # SPI MODEL
    rho_1_spi = CP.PropsSI('D', 'H', h_1, 'P', P_1, 'N2O') 
    m_dot_spi = Cd_hem_spi_dyer * A_inj_ox * np.sqrt( 2 * rho_1_spi * (P_1 - P_2)  )

    # HEM MODEL
    m_dot_hem = None
    downstream_pres_arr = np.linspace(P_2, P_1, 100)
    m_dot_hem_arr = []

    for pres in downstream_pres_arr:
        s_2 = CP.PropsSI('S', 'H', h_1, 'P', P_1, 'N2O') #assuming isentropic, upstream entropy equals downstream entropy
        h_2_hem = CP.PropsSI('H', 'S', s_2, 'P', pres, 'N2O')
        rho_2_hem = CP.PropsSI('D', 'S', s_2, 'P', pres, 'N2O')
        
        m_dot_hem = Cd_hem_spi_dyer * A_inj_ox * rho_2_hem * np.sqrt( 2 * np.abs(h_1 -  h_2_hem) )
        
        m_dot_hem_arr.append(m_dot_hem)

    m_dot_hem_crit = np.max(m_dot_hem_arr)
    P_crit = downstream_pres_arr[np.argmax(m_dot_hem_arr)]


    if P_2 < P_crit:
        #print("HEM predict choked flow: ", P_crit, m_dot_hem_crit)
        m_dot_hem = m_dot_hem_crit
        #P_2 = P_crit
        #plt.axvline( x = (P_1-P_2), color = 'rebeccapurple', linestyle = '-' )

    else:
        #print("HEM predict unchoked flow")
        
        s_1 = CP.PropsSI('S', 'H', h_1, 'P', P_1, 'N2O')
        h_2_hem = CP.PropsSI('H', 'S', s_1, 'P', P_2, 'N2O')
        rho_2_hem = CP.PropsSI('D', 'S', s_1, 'P', P_2, 'N2O')
        m_dot_hem = Cd_hem_spi_dyer * A_inj_ox * rho_2_hem * np.sqrt( 2 * (h_1 -  h_2_hem) )
        
        #plt.axvline( x = (P_1-P_2), color = 'honeydew', linestyle = '-' )
                                        

    # Dyer model

    #print(P_sat, CP.PropsSI('P', 'Q', 0, 'T', 282, 'N2O') )
    dyer_k = np.sqrt( (P_1 - P_2) / (P_sat - P_2) ) 
    #print(dyer_k, P_1, P_sat, P_2, P_crit)

    #NOTE: FOR THIS CASE REPLACED DYER CALL WITH SPI MODEL SINCE THAT WILL BE DOMINATE MASS FLOW PREDICITON
    #if(P_sat < P_2) and (m_dot_hem == m_dot_hem_crit) :
    #    print("(-) ", P_1, P_sat, P_2)


    m_dot_dyer = ((dyer_k/(1+dyer_k)) * m_dot_spi) + ((1/(1+dyer_k)) * m_dot_hem)
    
    #print(P_sat, P_2, m_dot_dyer, m_dot_hem, m_dot_spi)

    return m_dot_dyer







def proposed_model_inst(P_1, P_2, T_1): #TODO: ADD below constants TO MODEL INPUTS:
    
    #print(f"\n\n\nReduced Temp: {T_1/309.52}, Reduced Pres: {P_1/7.245e6}")
    
    
    all_err = 0.01

    #NOTE: GUESSING Cd
    Cd_hem_spi_dyer =  0.66
    Cd_high_supercharge = None #solving for
    A_inj_ox = 0.25*np.pi*((1.5e-3)**2) #m^2 #waxman test case

    # var declaration
    rho_1 = None
    h_1 = None
    x_1 = None
    G_sat = None
    m_dot = None

    # Start --> solve for two phase case at inlet

    #Honestly not sure what to do at this second, going to try
    # x = 0 since validating against subcooled test case but 
    # need to solve this for low subcooled state
    P_sat = CP.PropsSI('P', 'T', T_1, 'Q', 0, 'N2O')

    v_1_g = 1/CP.PropsSI('D', 'Q', 1, 'P', P_1, 'N2O')
    v_1_l = 1/CP.PropsSI('D', 'Q', 0, 'P', P_1, 'N2O')
    v_1_lg = v_1_g - v_1_l

    h_1_g = CP.PropsSI('H', 'Q', 1, 'P', P_1, 'N2O')
    h_1_l = CP.PropsSI('H', 'Q', 0, 'P', P_1, 'N2O')
    h_1_lg = h_1_g - h_1_l

    c_1_l = CP.PropsSI('CPMASS', 'Q', 0, 'P', P_1, 'N2O')

    
    rho_1 = CP.PropsSI('D', 'P', P_1, 'T', T_1, 'N2O')
    v_1 = 1/rho_1
    h_1 = CP.PropsSI('H', 'P', P_1, 'T', T_1, 'N2O')

    x_1 = 0 #ASSUMPTION, TESTING THIS, MIGHT NOT WORK, THATS OK


    omega_sat = (x_1*v_1_lg/v_1) + (c_1_l*T_1*P_1/v_1)*((v_1_lg/h_1_lg)**2) #i think this makes sense for it to be P_sat based off of Emerson's thesis but not sure anymore
    eta_sat = P_sat / P_1 
    

    #solve critical pressure
    #implicitly solve eta_crit
    
    eta_crit_sat = eta_sat #initial guess for critical pressure ratio
    while np.abs(TWOPHASEerror(eta_crit_sat, omega_sat) ) > all_err:
        eta_crit_sat = secant((lambda T: TWOPHASEerror(T, omega_sat)), eta_crit_sat)
    

    #print("eta sat: ", eta_sat, eta_crit_sat)
    P_crit_sat = eta_crit_sat * P_1
    #print("\n", eta_crit_sat * P_1, eta_crit_sat * P_sat, "\n")

    if P_crit_sat > P_2: 
        #print("SATURATED INLET CHOKED: ")
        G_sat = eta_crit_sat * np.sqrt(P_1*rho_1/omega_sat)

    else:
        #print("SATURATED INLET ***NOT*** CHOKED: ")
        G_sat = np.sqrt(P_1*rho_1) *np.sqrt( -2*(omega_sat * np.log(eta_sat) + (omega_sat -1)*(1 - eta_sat)) ) / (omega_sat*((1/eta_sat) - 1) + 1)


    #1) --> check inlet fluid state (sat liq. or subcooled)
    phase = CP.PropsSI('Phase', 'P', P_1, 'T', T_1, 'N2O')

    if phase == 5: #saturated liquid vapor
        #TODO: implement
        m_dot = None
    
    if phase == 0: #subcooled fluid 

        #setup subcooled fluid inlet 
        
        P_sat = CP.PropsSI('P', 'T', T_1, 'Q', 0, 'N2O')

        v_1_g = 1/CP.PropsSI('D', 'Q', 1, 'P', P_1, 'N2O')
        rho_1_l = CP.PropsSI('D', 'Q', 0, 'P', P_1, 'N2O')
        v_1_l = 1/rho_1_l
        v_1_lg = v_1_g - v_1_l

        h_1_g = CP.PropsSI('H', 'Q', 1, 'P', P_1, 'N2O')
        h_1_l = CP.PropsSI('H', 'Q', 0, 'P', P_1, 'N2O')
        h_1_lg = h_1_g - h_1_l

        c_1_l = CP.PropsSI('CPMASS', 'Q', 0, 'P', P_1, 'N2O') #BUG: ? assuming specific heat capacity at constant volume, thesis wasnt clear, might be a mistake

        rho_1 = CP.PropsSI('D', 'T', T_1, 'P', P_1, 'N2O')
        v_1 = 1/rho_1
        h_1 = CP.PropsSI('H', 'T', T_1, 'P', P_1, 'N2O')

        

        omega_sat = (c_1_l*T_1*P_sat/v_1_l)*( (v_1_lg/h_1_lg)**2)
        eta_transition =  2*omega_sat / (1 + 2*omega_sat)
        

        ### High subcooled
        if P_sat <= (eta_transition * P_1):
            #print("HIGH SUBCOOLED")

            P_sat = CP.PropsSI('P', 'T', T_1, 'Q', 0, 'N2O')
            #print("correction factor of 0.9")
            omega_sat = (c_1_l*T_1*P_sat/v_1_l)*( (v_1_lg/h_1_lg)**2)
            eta_transition =  2*omega_sat / (1 + 2*omega_sat)

            eta_crit = (P_sat / (P_1)) 
            P_crit = P_sat

            if P_2 < P_crit: #choked flow
                Cd_high_supercharge = 0.73
                m_dot = (Cd_high_supercharge*A_inj_ox) * np.sqrt( 2*(1-eta_crit)*P_1*rho_1)


                
            else: 
                #m_dot = dyer_model( Cd_hem_spi_dyer, A_inj_ox, P_1, P_sat, P_2, h_1 )
                
                #BUG? Emerson's thesis directs using dyer model here as ^, but the issue I have is that we are getting a negative denom in the dyer const K. Obv it is supposed to weight SPI model really high in this case so just using here
                #that might be a patchwork fix that covers something inherently wrong w the model
                # SPI MODEL
                rho_1_spi = CP.PropsSI('D', 'H', h_1, 'P', P_1, 'N2O') 
                m_dot = 1.105*Cd_hem_spi_dyer * A_inj_ox * np.sqrt( 2 * rho_1_spi * (P_1 - P_2)  )

                m_dot = dyer_model( Cd_hem_spi_dyer, A_inj_ox, P_1, P_sat, P_2, h_1 )
                #print("correction/smoothing factor of 1.105")
                
                #plt.axvline( x = (P_1-P_2), color = 'r', linestyle = '-' )
                print(P_1, P_sat, P_2)



        # Low subcooled
        else:

            ##print("low subcooled: ")

            ###NOTE: this seemed to fix low subcooled choked flow case
            eta = P_2 / P_1
            eta_crit_sat = eta #initial guess for critical pressure ratio

            while np.abs(LOWSUBCOOLEDerror(eta_crit_sat, eta_sat, omega_sat) ) > all_err:
                eta_crit_sat = secant((lambda T: LOWSUBCOOLEDerror(T, eta_sat, omega_sat)), eta_crit_sat)

            P_crit_low = eta_crit_sat * P_1


            if P_crit_low > P_2: #I had this equality mixed up for like a week fml... you remember when you are in elementary school, and they are teaching you about this
                #in terms of crocodile mouths going right or left, yea i clearly dont either send me back bruh
                #print("SATURATED INLET CHOKED: ")
                
                G_sat = eta_crit_sat * np.sqrt(P_1*rho_1/omega_sat)
                #print("saturated inlet choking: ", G_sat)
            else:
                #print("SATURATED INLET ***NOT*** CHOKED: ")
                G_sat = np.sqrt(P_1*rho_1) *np.sqrt( -2*(omega_sat * np.log(eta_sat) + (omega_sat -1)*(1 - eta_sat)) ) / (omega_sat*((1/eta_sat) - 1) + 1)

            #implicitly solve eta_crit_low
            eta = P_2 / P_1
            eta_crit_low = eta #initial guess for critical pressure ratio

            while np.abs(LOWSUBCOOLEDerror(eta_crit_low, eta_sat, omega_sat) ) > all_err:
                eta_crit_low = secant((lambda T: LOWSUBCOOLEDerror(T, eta_sat, omega_sat)), eta_crit_low)
            P_crit_low = eta_crit_low * P_1

            P_crit = eta_sat*P_crit_sat + (1-eta_sat)*P_crit_low #NOTE: using modified omega model to predict choking pressure and mass flow rate

            
            #check for choking
            if P_2 <= P_crit: #choked flow
                G_low = np.sqrt(rho_1_l * P_1) * np.sqrt( 2*(1-eta_sat) + 2*(omega_sat*eta_sat*np.log(eta_sat/eta_crit_low) - (omega_sat-1)*(eta_sat-eta_crit_low))) / (omega_sat*((eta_sat/eta_crit_low) - 1) + 1)
                m_dot= A_inj_ox *( (P_sat/P_1)*G_sat + (1-(P_sat/P_1))*G_low )

            else: #not choked use dyer model

                #solve choked: 
                G_sat_choked = eta_crit_sat * np.sqrt(P_1*rho_1/omega_sat)
                #implicitly solve eta_crit_low
                eta = P_2 / P_1
                eta_crit_low = eta #initial guess for critical pressure ratio

                while np.abs(LOWSUBCOOLEDerror(eta_crit_low, eta_sat, omega_sat) ) > all_err:
                    eta_crit_low = secant((lambda T: LOWSUBCOOLEDerror(T, eta_sat, omega_sat)), eta_crit_low)
                P_crit_low = eta_crit_low * P_1

                P_crit = eta_sat*P_crit_sat + (1-eta_sat)*P_crit_low #NOTE: using modified omega model to predict choking pressure and mass flow rate

                G_low_choked = np.sqrt(rho_1_l * P_1) * np.sqrt( 2*(1-eta_sat) + 2*(omega_sat*eta_sat*np.log(eta_sat/eta_crit_low) - (omega_sat-1)*(eta_sat-eta_crit_low))) / (omega_sat*((eta_sat/eta_crit_low) - 1) + 1)
                m_dot_choked = A_inj_ox *( (P_sat/P_1)*G_sat_choked + (1-(P_sat/P_1))*G_low_choked )

                #smoothing_factor = (6.083086*m_dot_choked + 0.844554896)
                #print("smoothing factor: ", smoothing_factor)
                #m_dot = (smoothing_factor)*dyer_model( Cd_hem_spi_dyer, A_inj_ox, P_1, P_sat, P_2, h_1 )
                m_dot = dyer_model( Cd_hem_spi_dyer, A_inj_ox, P_1, P_sat, P_2, h_1 )
                print(P_1, P_sat, P_2)

    return(m_dot)
    




#waxman_data = [ (280, 0.28e6, 'black'), (281, 0.55e6, 'red'), (282, 0.79e6, 'maroon' ), (281, 1.17e6, 'orangered'), (282, 1.42e6, 'pink'), (281, 1.85e6, 'yellow'), (282, 2e6, 'orange'), (283, 2.27e6, 'green'), (283, 2.55e6, 'blue') ]
#waxman_data = [(282, 0.79e6, 'maroon' ), (281, 1.17e6, 'orangered'), (282, 1.42e6, 'pink'), (281, 1.85e6, 'yellow'), (282, 2e6, 'orange'), (283, 2.27e6, 'green'), (283, 2.55e6, 'blue') ]
#waxman_data = [ (282, 0.79e6, 'maroon' ), (282, 2e6, 'orange'), ]
#waxman_data = [ (280, 0.28e6, 'black') , (281, 0.55e6, 'red')]
#waxman_data = [(282, 1.42e6, 'pink') ]
waxman_data = [ (283, 2.55e6, 'blue') ]

temp, P_super, correlating_colors = zip(*waxman_data)

for i in range(len(waxman_data)):
    T = temp[i]
    vapor_pressure = CP.PropsSI('P', 'T', T, 'Q', 0, 'N2O')  
    total_pres = vapor_pressure + P_super[i]

    P_upstream = total_pres
    P_downstream = 1e6
    P_arr = np.linspace(P_upstream-1e4 , P_downstream, 100)
    m_dot_arr = []
    delta_P_arr = []

    print(correlating_colors[i],"\n----------")

    for j in P_arr:
        x = proposed_model_inst(P_upstream, j,  T)
        m_dot_arr.append(x)
        delta_P_arr.append(P_upstream-j)

    print("\n\n\n")

    plt.plot(delta_P_arr,m_dot_arr, label = f'Supercharged {P_super[i]}', color = correlating_colors[i])

print("NOTE: correction factor of 0.9 in high subcooled P_sat case")
#plt.plot(delta_P_exp, m_dot_exp, label = 'exp data', color = 'pink')
plt.xlabel('delta P (MPa)')
plt.ylabel('Mass Flow Rate (kg/s)')
plt.title('delta P (MPa) vs Mass Flow Rate (kg/s)')
plt.grid(True)
plt.show()