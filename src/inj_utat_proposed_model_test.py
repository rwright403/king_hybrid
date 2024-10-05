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
    
    #print("function diff: ", function_diff, eta_crit)
    
    return function_diff

def LOWSUBCOOLEDerror(eta_crit, eta_sat, omega_sat):
    function_diff = (((omega_sat+(1/omega_sat)-2)/(2*eta_sat))*eta_crit**2) - (2*(omega_sat-1)*eta_crit) + (omega_sat*eta_sat*np.log(eta_crit/eta_sat)) + ((3/2)*omega_sat*eta_sat) - 1
    return function_diff

def proposed_model_inst(P_1, P_2, T_1):
    #constants and var declaration:
    all_err = 0.01

    #NOTE: GUESSING Cd
    Cd_ox =  0.6
    A_inj_ox = 0.25*np.pi*((1.5e-3)**2) #m^2

    rho_1 = None
    h_1 = None
    x_1 = None
    omega_sat = None
    G_sat = None
    m_dot = None


    # Start --> solve for two phase case at inlet

    #Honestly not sure what to do at this second, going to try
    # x = 1 since validating against subcooled test case but 
    # need to solve this for low subcooled state
    P_sat = CP.PropsSI('P', 'T', T_1, 'Q', 0, 'N2O')

    v_1_g = 1/CP.PropsSI('D', 'Q', 1, 'P', P_sat, 'N2O')
    v_1_l = 1/CP.PropsSI('D', 'Q', 0, 'P', P_sat, 'N2O')
    v_1_lg = v_1_g - v_1_l
    h_1_g = CP.PropsSI('H', 'Q', 1, 'P', P_sat, 'N2O')
    h_1_l = CP.PropsSI('H', 'Q', 0, 'P', P_sat, 'N2O')
    h_1_lg = h_1_g - h_1_l

    c_1_l = CP.PropsSI('CVMASS', 'Q', 0, 'P', P_sat, 'N2O') #BUG: ? assuming specific heat capacity at constant volume, thesis wasnt clear, might be a mistake

    x_1 = 0 #ASSUMPTION, TESTING THIS, MIGHT NOT WORK, THATS OK

    rho_1 = CP.PropsSI('D', 'Q', x_1, 'T', T_1, 'N2O')
    v_1 = 1/rho_1
    h_1 = CP.PropsSI('H', 'Q', x_1, 'T', T_1, 'N2O')

    omega_sat = (x_1*v_1_lg/v_1) + (c_1_l*T_1*P_sat/v_1)*((v_1_lg/h_1_lg)**2)
    eta_sat = P_2 / P_sat 

    #solve critical pressure
    #implicitly solve eta_crit
    eta_crit_sat = eta_sat #initial guess for critical pressure ratio
    while np.abs(TWOPHASEerror(eta_crit_sat, omega_sat) ) > all_err:
        eta_crit_sat = secant((lambda T: TWOPHASEerror(T, omega_sat)), eta_crit_sat)
    
    P_crit_sat = eta_crit_sat * P_sat

    if P_crit_sat < P_2:
        G_sat = eta_crit_sat * np.sqrt(P_1*rho_1/omega_sat)
    else:
        G_sat = np.sqrt( -2*(omega_sat * np.log(eta_sat) + (omega_sat -1)*(1 - eta_sat)) ) / (omega_sat*((1/eta_sat) - 1) + 1)
    


    #1) --> check inlet fluid state (sat liq. or subcooled)
    phase = CP.PropsSI('Phase', 'P', P_1, 'T', T_1, 'N2O')

    if phase == 5: #saturated
        #NOTE: DOES THIS PREDICT CHOKING OR NOT?
        m_dot = Cd_ox*A_inj_ox * G_sat
    
    if phase == 0: #subcooled fluid 

        #setup subcooled fluid inlet 
        #NOTE: omega model EOS assumes initial pressure is P_sat corresponding to T_1, even for subcooled conditions
        P_sat = CP.PropsSI('P', 'T', T_1, 'Q', 0, 'N2O')

        v_1_g = 1/CP.PropsSI('D', 'Q', 1, 'P', P_sat, 'N2O')
        v_1_l = 1/CP.PropsSI('D', 'Q', 0, 'P', P_sat, 'N2O')
        v_1_lg = v_1_g - v_1_l

        h_1_g = CP.PropsSI('H', 'Q', 1, 'P', P_sat, 'N2O')
        h_1_l = CP.PropsSI('H', 'Q', 0, 'P', P_sat, 'N2O')
        h_1_lg = h_1_g - h_1_l

        c_1_l = CP.PropsSI('CVMASS', 'Q', 0, 'P', P_sat, 'N2O') #BUG: ? assuming specific heat capacity at constant volume, thesis wasnt clear, might be a mistake

        rho_1 = CP.PropsSI('D', 'P', P_1, 'T', T_1, 'N2O')
        v_1 = 1/rho_1
        h_1 = CP.PropsSI('H', 'P', P_1, 'T', T_1, 'N2O')

        omega_sat = (c_1_l*T_1*P_sat/v_1)*((v_1_lg/h_1_lg)**2)
        eta_transition =  2 * omega_sat / ( 1 + 2*omega_sat)


        #eta_sat = P_sat / P_1 NOTE: ^ AND THIS DONT CORRELATE
        #print("checking someth: ", eta_sat,  2 * omega_sat / ( 1 + 2*omega_sat))
        

        m_dot = None

        #is assuming flow will always be choked here wrong?

        # High subcooled
        if P_sat < (eta_transition * P_1):

            #at high subcooling conditions, no vapor formed during flow through nozzle:
            # crit pressure ratio = P_sat / P_inlet
            #P_crit_high = P_sat
            eta_crit = P_sat / P_1

            #i think for modified omega if its high supercharged we can assume its choking according to Emerson's thesis
            m_dot = (Cd_ox*A_inj_ox) * np.sqrt( 2*(1-eta_crit)*P_1*rho_1)

            """
            #check for choking
            print(P_2, P_crit_high)
            if P_2 <= P_crit_high: #choked flow
                #m_dot = (Cd_ox*A_inj_ox) * np.sqrt( 2*(1-eta_crit)*P_1*rho_1)
                m_dot = (Cd_ox*A_inj_ox) * rho_1 * np.sqrt(2*P_1 * (P_1 - P_crit_high))
            else:
                eta = P_2 / P_1
                #NOTE: PROBLEM WITH UNCHOKED HIGH SUBCOOLED FORMULA
                #m_dot = (Cd_ox*A_inj_ox) * np.sqrt( 2*(1-eta_sat) + 2*(omega_sat*eta_sat*np.log(eta_sat/eta) - (omega_sat-1)*(eta_sat-eta))) / (omega_sat*((eta_sat/eta) - 1) + 1)
                #testing someth

                #for chamber pressures above the choking pressure, the dyer model is used and connected to the choked regime using a blending function

                #solving "modified omega crit pressure" then using that to solve choked flow
            """
            print("high supercharge m_dot_ox: ", m_dot)


        # Low subcooled
        if P_sat > (eta_transition * P_1):
            print("low subcooled: ", P_1)

            #NOTE: using modified omega model to predict choking pressure and mass flow rate
            #for all p_cc below choking (for all choked flow) use

            #implicitly solve eta_crit_low
            eta = P_2 / P_1
            eta_crit_low = eta #initial guess for critical pressure ratio

            while np.abs(LOWSUBCOOLEDerror(eta_crit_low, eta_sat, omega_sat) ) > all_err:
                eta_crit_low = secant((lambda T: LOWSUBCOOLEDerror(T, eta_sat, omega_sat)), eta_crit_low)

            #NOTE: This is smoothed
            P_crit_low = eta_crit_low * P_1

            #this is the modified omega model
            P_crit = eta_sat*P_crit_sat + (1-eta_sat)*P_crit_low
            #print(P_crit_sat, P_crit_low)

        
            #check for choking
            if P_2 <= P_crit: #choked flow



                G_low = np.sqrt( 2*(1-eta_sat) + 2*(omega_sat*eta_sat*np.log(eta_sat/eta_crit_low) - (omega_sat-1)*(eta_sat-eta_crit_low))) / (omega_sat*((eta_sat/eta_crit_low) - 1) + 1)
                m_dot = (Cd_ox*A_inj_ox) * (eta_sat * G_sat + (1-eta_sat)*G_low)
                #print("low subcooled choked: ", m_dot, eta_sat * G_sat, (1-eta_sat)*G_low )

            else: #solve dyer model and connect with blending function
                G_low = np.sqrt( 2*(1-eta_sat) + 2*(omega_sat*eta_sat*np.log(eta_sat/eta) - (omega_sat-1)*(eta_sat-eta))) / (omega_sat*((eta_sat/eta) - 1) + 1)
                #using dyer model here?

                #DYER MODEL
                # SPI MODEL
                rho_2_spi = CP.PropsSI('D', 'H', h_1, 'P', P_1, 'N2O') #is isentropic valid for this model?
                m_dot_spi = Cd_ox * A_inj_ox * np.sqrt( 2 * rho_2_spi * (P_1 - P_2)  )

                # HEM MODEL
                s_2 = CP.PropsSI('S', 'H', h_1, 'P', P_1, 'N2O') #assuming isentropic, upstream entropy equals downstream entropy
                h_2_hem = CP.PropsSI('H', 'S', s_2, 'P', P_2, 'N2O')
                rho_2_hem = CP.PropsSI('D', 'S', s_2, 'P', P_2, 'N2O')

                m_dot_hem = Cd_ox * A_inj_ox * rho_2_hem * np.sqrt( 2 * (h_1 -  h_2_hem) )
                                
                # Dyer MODEL 
                #print("low subcooled not choked: dyer k denom: ", (P_sat - P_2) )
                dyer_k = 0 #np.sqrt( (P_1 - P_2) / (P_sat - P_2) ) #negative denomenator
                #if P_2 > P_sat... then it would be exiting as a subcooled liquid  
                # doesnt this violate when hem would be used?  
                 
                m_dot = ((dyer_k/(1+dyer_k)) * m_dot_spi) + ((1/(1+dyer_k)) * m_dot_hem)
                #print(m_dot, P_sat ,(eta_transition * P_1))

                #m_dot = (P_sat/P_1)*((Cd_ox*A_inj_ox)*G_sat)+ ((P_sat/P_1))*m_dot_dyer #THIS IS WRONG!
                
                #print("low subcooled not choked: ", m_dot)
                ##print("low subcooled: ", m_dot)

    #print(m_dot)
    return(m_dot)



P_arr = np.linspace(4.36e6, 6.93e6, 100)
m_dot_arr = []

for i in P_arr:
    if i != 4.36e6:
        x = proposed_model_inst(i, 4.36e6, 282)
        #print(x)
        m_dot_arr.append(x)
    else:
        m_dot_arr.append(0)


"""
P_arr = np.linspace( 6.8e6, 4.36e6, 100)
m_dot_arr = []

for i in P_arr:
    if i != 4.36e6:
        print ("delta p: ", 6.8e6-i)
        x = proposed_model_inst(6.93e6, i, 282)
        #print(x)
        m_dot_arr.append(x)
    else:
        m_dot_arr.append(0)
"""

delta_P_arr = []
for i in P_arr:
    delta_P_arr.append(i-4.36e6)

plt.plot(delta_P_arr,m_dot_arr)
plt.xlabel('delta P (MPa)')
plt.ylabel('Mass Flow Rate (kg/s)')
plt.title('delta P (MPa) vs Mass Flow Rate (kg/s)')
plt.grid(True)
plt.show()

