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
#P_REF = 401325 #Pa

n2o = Chemical('nitrous oxide', T=T_REF)

MW = (n2o.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o.omega - 0.26992*n2o.omega**2

TANK_DIAM = 0.0254*5.5 #m
CS_AREA = 0.25*np.pi*(TANK_DIAM**2) #m^2
g = 9.81 #m/s^2

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

"""
def hem_predict_choking(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, h_1):
    m_dot_hem = None
    downstream_pres_arr = np.linspace(P_2, P_1, 100)
    m_dot_hem_arr = []

    for pres in downstream_pres_arr:
        s_2 = CP.PropsSI('S', 'H', h_1, 'P', P_1, 'N2O') #assuming isentropic, upstream entropy equals downstream entropy
        h_2_hem = CP.PropsSI('H', 'S', s_2, 'P', pres, 'N2O')
        rho_2_hem = CP.PropsSI('D', 'S', s_2, 'P', pres, 'N2O')
            
        m_dot_hem = Cd_hem_spi_dyer * A_inj_ox * rho_2_hem * np.sqrt( 2 * np.abs(h_1 -  h_2_hem) )
            
        m_dot_hem_arr.append(m_dot_hem)

    P_crit = downstream_pres_arr[np.argmax(m_dot_hem_arr)]

    #print(P_2, P_crit)

    if P_2 < P_crit: #flow is choked
        return P_crit
    else: #flow is unchoked
        return P_2
"""

def spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit, h_1):
    m_dot_spi = Cd_hem_spi_dyer * A_inj_ox * np.sqrt( 2 * rho_tank_exit * (P_1 - P_2)  )
    print(rho_tank_exit)
    return m_dot_spi

def hem_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, h_1):
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

    if P_2 < P_crit: #flow is choked
        m_dot_hem = m_dot_hem_crit
        print("HEM PREDICT CHOKING!!!! ITS CHOKING!!!!")

    else: #flow is unchoked
        s_1 = CP.PropsSI('S', 'H', h_1, 'P', P_1, 'N2O')
        h_2_hem = CP.PropsSI('H', 'S', s_1, 'P', P_2, 'N2O')
        rho_2_hem = CP.PropsSI('D', 'S', s_1, 'P', P_2, 'N2O')
        m_dot_hem = Cd_hem_spi_dyer * A_inj_ox * rho_2_hem * np.sqrt( 2 * (h_1 -  h_2_hem) )

    return m_dot_hem


def dyer_model( Cd_hem_spi_dyer, A_inj_ox, P_1, P_sat, P_2, rho_tank_exit, h_1):

    m_dot_spi = spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit, h_1)
                                        
    #NOTE: FOR THIS CASE NO CAVITATION SO WE ARE JUST USING THE SPI MODEL
    if(P_sat < P_2):
        m_dot_dyer = m_dot_spi

    #NOTE: ELSE TWO PHASE AT INJ OUTLET AND USE DYER TO ACCOUNT FOR TWO PHASE EFFECTS
    else:
        m_dot_hem = hem_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, h_1)
        print("HEM: ", m_dot_hem)
        
        dyer_k = np.sqrt( (P_1 - P_2) / (P_sat - P_2) ) 
        m_dot_dyer = ((dyer_k/(1+dyer_k)) * m_dot_spi) + ((1/(1+dyer_k)) * m_dot_hem)
    
    return m_dot_dyer


#NOTE: not using subcooled inlet stuff in this tank model
def modified_emerson_and_mohammad_model_inst(P_1, P_2, T_1, x_1, A_inj_1, rho_1, h_1): #TODO: ADD below constants TO MODEL INPUTS:
    all_err = 0.005 #TODO: maybe UPDATE TO USE WHAT THE REST OF THE FUNCTION USES?



    #1) --> check inlet fluid state (sat liq. or subcooled)
    phase = CP.PropsSI('Phase', 'P', P_1, 'D', rho_1, 'N2O')

    #print(phase, CP.PhaseSI('P', P_1, 'D', rho_1, 'N2O'))



    #NOTE: IF SUBCOOLED LIQUID x_1 SHOULD BE 0 ASSUMPTION, TESTING THIS, MIGHT NOT WORK, THATS OK TODO: MAKE AN INPUT to this function

    Cd_hem_spi_dyer =  0.66 #EDUCATED GUESS FOR CD BASED ON SPREADSHEET, COULD IMPLEMENT THOSE EQNS TO IMPROVE ACCURACY BUT THIS IS LIKELY OK
    Cd_high_supercharge = None #solving for

    # var declaration
    G_sat = None
    m_dot = None

    # Start --> solve for two phase case at inlet

    P_sat = CP.PropsSI('P', 'T', T_1, 'Q', 0, 'N2O')

    v_1_g = 1/CP.PropsSI('D', 'Q', 1, 'P', P_1, 'N2O')
    v_1_l = 1/CP.PropsSI('D', 'Q', 0, 'P', P_1, 'N2O')
    v_1_lg = v_1_g - v_1_l
    h_1_g = CP.PropsSI('H', 'Q', 1, 'P', P_1, 'N2O')
    h_1_l = CP.PropsSI('H', 'Q', 0, 'P', P_1, 'N2O')
    h_1_lg = h_1_g - h_1_l

    c_1_l = CP.PropsSI('CPMASS', 'Q', 0, 'P', P_1, 'N2O')

    #rho_1 = CP.PropsSI('D', 'P', P_1, 'T', T_1, 'N2O')
    v_1 = 1/rho_1
    #h_1 = CP.PropsSI('H', 'P', P_1, 'T', T_1, 'N2O')

    #1) --> check inlet fluid state (sat liq. or subcooled)
    #phase = CP.PropsSI('Phase', 'P', P_1, 'T', T_1, 'N2O')

    omega_sat = (x_1*v_1_lg/v_1) + (c_1_l*T_1*P_sat/v_1)*((v_1_lg/h_1_lg)**2) #i think this makes sense for it to be P_sat based off of Emerson's thesis but not sure anymore
    eta_sat = P_sat / P_1 

    #solve critical pressure
    #implicitly solve eta_crit
    
    eta_crit_sat = eta_sat #initial guess for critical pressure ratio
    while np.abs(TWOPHASEerror(eta_crit_sat, omega_sat) ) > all_err:
        eta_crit_sat = secant((lambda T: TWOPHASEerror(T, omega_sat)), eta_crit_sat)
    

    #print("eta sat: ", eta_sat, eta_crit_sat)
    P_crit_sat = eta_crit_sat * P_1
    #print("\n", eta_crit_sat * P_1, eta_crit_sat * P_sat, "\n")

    if P_crit_sat >= P_2: #Saturated Inlet Choked
        G_sat = eta_crit_sat * np.sqrt(P_1*rho_1/omega_sat)

        print("sat inlet choked")

    else: #Saturated Inlet Not Choked
        pratio = P_2/P_sat
        G_sat = np.sqrt(P_1*rho_1) *np.sqrt( -2*(omega_sat * np.log(pratio) + (omega_sat -1)*(1 - pratio)) ) / (omega_sat*((1/(pratio)) - 1) + 1)



        #print("sat inlet not choked", G_sat,  (-2*(omega_sat * np.log(eta_sat))) , ((omega_sat -1)*(1 - eta_sat)), omega_sat, eta_sat)

    if phase == 6: #saturated liquid vapor
        k_cavitation_const = (P_1 - P_sat) / (P_1 - P_2)
        if k_cavitation_const == 0:
            k_cavitation_const = 1 #cavitation constant only valid for P_sat > P_2, otherwise K = 1
        Cd_twophase = 0.386 + 0.361*np.sqrt(k_cavitation_const)
        m_dot =  ( Cd_twophase *A_inj_1) * G_sat  

        
        #print("line 232", m_dot, rho_1, P_1, eta_crit_sat, omega_sat, x_1)
        #checked A_inj_1 is correct, fixed issue with cavitation constant k

    
    if phase == 0: #subcooled fluid 

        #setup subcooled fluid inlet 
        
        rho_1_l = CP.PropsSI('D', 'Q', 0, 'P', P_1, 'N2O')
        eta_transition =  2*omega_sat / (1 + 2*omega_sat)

        ### High subcooled
        if P_sat <= (eta_transition * P_1):
            #print("HIGH SUBCOOLED")

            P_sat = 0.9*CP.PropsSI('P', 'T', T_1, 'Q', 0, 'N2O')
            #print("correction factor of 0.9")
            omega_sat = (c_1_l*T_1*P_sat/v_1_l)*( (v_1_lg/h_1_lg)**2) #TODO: put corr factor on omega instead?
            eta_transition =  2*omega_sat / (1 + 2*omega_sat)

            eta_crit = (P_sat / (P_1)) 
            P_crit = P_sat

            if P_2 < P_crit: #choked flow
                Cd_high_supercharge = 0.73
                m_dot = (Cd_high_supercharge*A_inj_1) * np.sqrt( 2*(1-eta_crit)*P_1*rho_1)

            else: 
                m_dot = 1.105*dyer_model( Cd_hem_spi_dyer, A_inj_1, P_1, P_sat, P_2, rho_1, h_1 )
                #print("correction/smoothing factor of 1.105")




        ### Low subcooled
        else:
            ###NOTE: this seemed to fix low subcooled choked flow case
            eta = P_2 / P_1
            eta_crit_sat = eta #initial guess for critical pressure ratio

            while np.abs(LOWSUBCOOLEDerror(eta_crit_sat, eta_sat, omega_sat) ) > all_err:
                eta_crit_sat = secant((lambda T: LOWSUBCOOLEDerror(T, eta_sat, omega_sat)), eta_crit_sat)

            P_crit_low = eta_crit_sat * P_1


            if P_crit_low > P_2: #saturated inlet choked
                G_sat = eta_crit_sat * np.sqrt(P_1*rho_1/omega_sat)

            else: #saturated inlet unchoked
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
                m_dot= A_inj_1 *( (P_sat/P_1)*G_sat + (1-(P_sat/P_1))*G_low )

            else: #not choked use dyer model --> smoothing factor based on choked mass flow rate!

                #solve choked: 
                G_sat_choked = eta_crit_sat * np.sqrt(P_1*rho_1/omega_sat)
                #implicitly solve eta_crit_low
                eta = P_2 / P_1
                eta_crit_low = eta #initial guess for critical pressure ratio

                while np.abs(LOWSUBCOOLEDerror(eta_crit_low, eta_sat, omega_sat) ) > all_err:
                    eta_crit_low = secant((lambda T: LOWSUBCOOLEDerror(T, eta_sat, omega_sat)), eta_crit_low)
                
                G_low_choked = np.sqrt(rho_1_l * P_1) * np.sqrt( 2*(1-eta_sat) + 2*(omega_sat*eta_sat*np.log(eta_sat/eta_crit_low) - (omega_sat-1)*(eta_sat-eta_crit_low))) / (omega_sat*((eta_sat/eta_crit_low) - 1) + 1)
                m_dot_choked = A_inj_1 *( (P_sat/P_1)*G_sat_choked + (1-(P_sat/P_1))*G_low_choked )
                
                smoothing_factor = (6.083086 * m_dot_choked + 0.8304554896) #I made this up to get the transition to be smooth
                #print("smoothing factor: ", smoothing_factor)
                m_dot = (smoothing_factor)*dyer_model( Cd_hem_spi_dyer, A_inj_1, P_1, P_sat, P_2, rho_1, h_1 )

    return(m_dot)


P_atm = 1e5 #Pa
T_atm = 286.5 #K
rho_atm = 1.225 #kg/m^3

m_nos = 20 #kg
P_tank = 45e5 #Pa
V_tank = 0.0354 #m^3

diam_out = 0.230 #m #NOTE: thesis didn't provide tank geometry, estimated based off of G type nos dimensions (approx equivalent mass to Karabeyoglu run tank)
diam_in = 0.215 #m
rho_wall = 2770 #kg/m^3
k_w = 237 #W/(m K)

Cd_1 = 0.425
A_inj_1 = 0.00003 #m^3 NOTE: GUESS
P_cc = 1.03e6 #Pa

#modified_emerson_and_mohammad_model_inst(P_tank, P_cc, T_atm, x_1, A_inj_1, rho_1, h_1)


P_1 = 4.5e6 #Pa
rho_1 = CP.PropsSI("D", "P", P_1, "Q", 0, "N2O") - 20 #kg/m^3
T_1 = CP.PropsSI("T", "P", P_1, "Q", 0, "N2O") #k



phase_test_1 = CP.PhaseSI('P', P_1, 'D', rho_1, 'N2O')
h_test_1 = CP.PropsSI('H', 'P', P_1, 'D', rho_1, 'N2O')
phase_test_2 = CP.PhaseSI('P', T_1, 'D', rho_1, 'N2O')
h_test_2 = CP.PropsSI('H', 'P', T_1, 'D', rho_1, 'N2O')


print("1 vs 2 phase: ", phase_test_1, phase_test_2)
print("1 vs 2  enthalpy: ", h_test_1, h_test_2)


n2o_ig = Chemical('N2O', T=T_1) 
preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_1, P=P_tank)
h_liq = preos_l.H_dep_l/MW + n2o_ig.Cpg*(T_1 - T_REF)

print("thermo enthalpy: ", h_liq)


#setup
def zilliac2005_h(T):
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
    return F_1 + val 