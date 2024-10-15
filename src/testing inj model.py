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



def dyer_model( Cd_hem_spi_dyer, A_inj_ox, P_1, P_sat, P_2, h_1):

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
    dyer_k = np.sqrt( (P_1 - P_2) / (P_sat - P_2) ) 
    print(dyer_k, P_1, P_sat, P_2, P_crit)

    #NOTE: FOR THIS CASE REPLACED DYER CALL WITH SPI MODEL SINCE THAT WILL BE DOMINATE MASS FLOW PREDICITON
    if(P_sat < P_2) and (m_dot_hem == m_dot_hem_crit) :
        print("(-) ", P_1, P_sat, P_2)


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

    c_1_l = CP.PropsSI('CPMASS', 'Q', 0, 'P', P_1, 'N2O') #BUG: ? assuming specific heat capacity at constant pressure, thesis wasnt clear, might be a mistake

    

    rho_1 = CP.PropsSI('D', 'P', P_1, 'T', T_1, 'N2O')
    v_1 = 1/rho_1
    h_1 = CP.PropsSI('H', 'P', P_1, 'T', T_1, 'N2O')

    x_1 = 0 #ASSUMPTION, TESTING THIS, MIGHT NOT WORK, THATS OK

    #NOTE: MIGHT NEED TO USE P_1 UP THERE ^ this was changed to help fix choked subcooled, negligible effect tho, if other errors introduced, try changing back but i dont think this is a likely cause

    #1) --> check inlet fluid state (sat liq. or subcooled)
    phase = CP.PropsSI('Phase', 'P', P_1, 'T', T_1, 'N2O')

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

    if P_crit_sat > P_2: #I had this equality mixed up for like a week fml... you remember when you are in elementary school, and they are teaching you about this
        #in terms of crocodile mouths going right or left, yea i clearly dont either send me back bruh
        #print("SATURATED INLET CHOKED: ")
        
        
        G_sat = eta_crit_sat * np.sqrt(P_1*rho_1/omega_sat)
        #print("saturated inlet choking: ", G_sat)
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
        #NOTE: CHANGED TO CPMASS SEEMS TO FIX 0.28 test case, NEED TO CHECK, REALLY NEED TO CHECK AND NOT ASSUME

        rho_1 = CP.PropsSI('D', 'T', T_1, 'P', P_1, 'N2O')
        v_1 = 1/rho_1
        h_1 = CP.PropsSI('H', 'T', T_1, 'P', P_1, 'N2O')

        #omega_sat here for an initially subcooled fluid (NO VAPOR TERM)
        omega_sat = (c_1_l*T_1*P_sat/v_1_l)*( (v_1_lg/h_1_lg)**2)
        eta_transition =  2*omega_sat / (1 + 2*omega_sat)
        

        ### High subcooled
        if P_sat <= (eta_transition * P_1):
            #print("HIGH SUBCOOLED")

            P_sat = 0.9*CP.PropsSI('P', 'T', T_1, 'Q', 0, 'N2O')
            print("additonal factor of 0.9 for debugging")
            #omega_sat here for an initially subcooled fluid (NO VAPOR TERM)
            omega_sat = (c_1_l*T_1*P_sat/v_1_l)*( (v_1_lg/h_1_lg)**2)
            eta_transition =  2*omega_sat / (1 + 2*omega_sat)

            eta_crit = (P_sat / (P_1)) #*0.89 #fixes everything
            P_crit = P_sat

            if P_2 < P_crit: #choked flow
                Cd_high_supercharge = 0.73
                m_dot = (Cd_high_supercharge*A_inj_ox) * np.sqrt( 2*(1-eta_crit)*P_1*rho_1)
                print("here!!! ", m_dot)
                #print( m_dot, eta_crit, P_1, (P_1-P_crit), rho_1, 1/v_1_l)
                #m_dot = (Cd_high_supercharge*A_inj_ox) * np.sqrt(2*rho_1*(P_1-P_sat))

                #print(eta_crit, eta_crit*0.89, eta_crit_sat)


                
            else: 
                #m_dot = dyer_model( Cd_hem_spi_dyer, A_inj_ox, P_1, P_sat, P_2, h_1 )
                
                #BUG? Emerson's thesis directs using dyer model here as ^, but the issue I have is that we are getting a negative denom in the dyer const K. Obv it is supposed to weight SPI model really high in this case so just using here
                #that might be a patchwork fix that covers something inherently wrong w the model
                # SPI MODEL
                rho_1_spi = CP.PropsSI('D', 'H', h_1, 'P', P_1, 'N2O') 
                m_dot = Cd_hem_spi_dyer * A_inj_ox * np.sqrt( 2 * rho_1_spi * (P_1 - P_2)  )
                
                #plt.axvline( x = (P_1-P_2), color = 'r', linestyle = '-' )



        # Low subcooled
        else:

            ##print("low subcooled: ")

            ###NOTE: this seemed to fix low subcooled choked flow case
            eta = P_2 / P_1
            #BUG????? ETA CRIT LOW HERE THEN SAT BELOW???? SEEMS TO WORK
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

            #this is the modified omega model
            P_crit = eta_sat*P_crit_sat + (1-eta_sat)*P_crit_low #NOTE: using modified omega model to predict choking pressure and mass flow rate

            
            #check for choking
            if P_2 <= P_crit: #choked flow
                G_low = np.sqrt(rho_1_l * P_1) * np.sqrt( 2*(1-eta_sat) + 2*(omega_sat*eta_sat*np.log(eta_sat/eta_crit_low) - (omega_sat-1)*(eta_sat-eta_crit_low))) / (omega_sat*((eta_sat/eta_crit_low) - 1) + 1)

                #print("checking G_low: ", G_low, np.sqrt(rho_1_l * P_1) * np.sqrt( 2*(1-eta_sat) + 2*(omega_sat*eta_sat*np.log(eta_sat/eta) - (omega_sat-1)*(eta_sat-eta))) / (omega_sat*((eta_sat/eta) - 1) + 1))

                #m_dot = A_inj_ox*G_low  #A_inj_ox *( (P_sat/P_1)*G_sat + (1-(P_sat/P_1))*G_low )
                m_dot = A_inj_ox *( (P_sat/P_1)*G_sat + (1-(P_sat/P_1))*G_low )

                #m_dot = None
                #print(G_sat, G_low)
                #print(eta_crit_low, eta_crit_sat)

            else: #not choked use dyer model
                m_dot = dyer_model( Cd_hem_spi_dyer, A_inj_ox, P_1, P_sat, P_2, h_1 )

    return(m_dot)
    




waxman_data = [ (280, 0.28e6, 'black'), (281, 0.55e6, 'red'), (282, 0.79e6, 'maroon' ), (281, 1.17e6, 'orangered'), (282, 1.42e6, 'pink'), (281, 1.85e6, 'yellow'), (282, 2e6, 'orange'), (283, 2.27e6, 'green'), (283, 2.55e6, 'blue') ]
#waxman_data = [(282, 0.79e6, 'maroon' ), (281, 1.17e6, 'orangered'), (282, 1.42e6, 'pink'), (281, 1.85e6, 'yellow'), (282, 2e6, 'orange'), (283, 2.27e6, 'green'), (283, 2.55e6, 'blue') ]
#waxman_data = [ (282, 0.79e6, 'maroon' ), (282, 2e6, 'orange'), ]
#waxman_data = [ (280, 0.28e6, 'black') , (281, 0.55e6, 'red')]
#waxman_data = [(282, 1.42e6, 'pink') ]
#waxman_data = [ (283, 2.55e6, 'blue') ]

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
#plt.plot(delta_P_exp, m_dot_exp, label = 'exp data', color = 'pink')
plt.xlabel('delta P (MPa)')
plt.ylabel('Mass Flow Rate (kg/s)')
plt.title('delta P (MPa) vs Mass Flow Rate (kg/s)')
plt.grid(True)
plt.show()