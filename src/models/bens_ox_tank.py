#Adapted from Benjamin Klammer's MuleSim3 Matlab Script --> shoutout

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

def Verror(T_tank, U_tank, m_ox, V_tank):

    rho_liq = CP.PropsSI('D', 'Q', 0, 'T', T_tank, 'N2O')  # Density of liquid nitrous (kg/m^3)
    rho_vap = CP.PropsSI('D', 'Q', 1, 'T', T_tank, 'N2O')  # Density of nitrous gas (kg/m^3)

    u_liq = CP.PropsSI('U', 'Q', 0, 'T', T_tank, 'N2O')  # Internal Energy of liquid nitrous (J/kg)
    u_vap = CP.PropsSI('U', 'Q', 1, 'T', T_tank, 'N2O')  # Internal Energy of nitrous gas (J/kg)

    x_tank = (U_tank/m_ox - u_liq)/(u_vap - u_liq)
    V_tank_diff = m_ox * ( ((1-x_tank) /rho_liq) + (x_tank/rho_vap) ) - V_tank
    return V_tank_diff

def uerror(T, rho_tank, u_tank): #TODO:(u_tank + 7.3397e+5)
    u2 = thermo_span_wagner(rho_tank, T, 'u') + 7.3397e+5
    u_diff = u2 - u_tank
    return u_diff

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

def TWOPHASEerror(eta_crit, omega):
    function_diff = (eta_crit**2) + ((omega**2)-2*omega)*((1-eta_crit)**2) + 2*(omega**2)*np.log(eta_crit) + 2*(omega**2)*(1-eta_crit)  
    return function_diff

def LOWSUBCOOLEDerror(eta_crit, eta_sat, omega_sat):
    function_diff = (((omega_sat+(1/omega_sat)-2)/(2*eta_sat))*(eta_crit**2)) - (2*(omega_sat-1)*eta_crit) + (omega_sat*eta_sat*np.log(eta_crit/eta_sat)) + ((3/2)*omega_sat*eta_sat) - 1
    return function_diff


def spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit): 
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

    else: #flow is unchoked
        s_1 = CP.PropsSI('S', 'H', h_1, 'P', P_1, 'N2O')
        h_2_hem = CP.PropsSI('H', 'S', s_1, 'P', P_2, 'N2O')
        rho_2_hem = CP.PropsSI('D', 'S', s_1, 'P', P_2, 'N2O')
        m_dot_hem = Cd_hem_spi_dyer * A_inj_ox * rho_2_hem * np.sqrt( 2 * (h_1 -  h_2_hem) )

    return m_dot_hem


def dyer_model( Cd_hem_spi_dyer, A_inj_ox, P_1, P_sat, P_2, rho_tank_exit, h_1):

    m_dot_spi = spi_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, rho_tank_exit)
                                        
    #NOTE: FOR THIS CASE NO CAVITATION SO WE ARE JUST USING THE SPI MODEL
    if(P_sat < P_2):
        m_dot_dyer = m_dot_spi

    #NOTE: ELSE TWO PHASE AT INJ OUTLET AND USE DYER TO ACCOUNT FOR TWO PHASE EFFECTS
    else:
        m_dot_hem = hem_model(Cd_hem_spi_dyer, A_inj_ox, P_1, P_2, h_1)
        
        dyer_k = np.sqrt( (P_1 - P_2) / (P_sat - P_2) ) 
        m_dot_dyer = ((dyer_k/(1+dyer_k)) * m_dot_spi) + ((1/(1+dyer_k)) * m_dot_hem)
    
    return m_dot_dyer


#NOTE: not using subcooled inlet stuff in this tank model
def modified_emerson_and_mohammad_model_inst(P_1, P_2, T_1, x_1, A_inj_1, rho_1, h_1): #TODO: ADD below constants TO MODEL INPUTS:
    all_err = 0.01 #TODO: maybe UPDATE TO USE WHAT THE REST OF THE FUNCTION USES?



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




class model():
    def __init__(self, oxidizer, timestep, m_ox, Cd_1, A_inj_1, V_tank, P_tank, P_cc, all_error, inj_model):
        self.oxidizer = oxidizer
        self.timestep = timestep
        self.m_ox = m_ox
        self.Cd_1 = Cd_1
        self.A_inj_1 = A_inj_1
        self.m_dot_ox = 0
        self.V_tank = V_tank
        self.P_tank = P_tank
        self.P_cc = P_cc
        self.all_error = all_error
        self.inj_model = inj_model




        #setup 
        self.t = 0
        self.m_dot_ox_prev = self.m_dot_ox

        #setup thermo data
        self.rho_liq = CP.PropsSI('D', 'Q', 0, 'P', self.P_tank, 'N2O')  # Density of liquid nitrous (kg/m^3)
        self.rho_vap = CP.PropsSI('D', 'Q', 1, 'P', self.P_tank, 'N2O')  # Density of nitrous gas (kg/m^3)

        self.u_liq = CP.PropsSI('U', 'Q', 0, 'P', self.P_tank, 'N2O')  # Internal Energy of liquid nitrous (kJ/kg)
        self.u_vap = CP.PropsSI('U', 'Q', 1, 'P', self.P_tank, 'N2O')  # Internal Energy of nitrous gas (kJ/kg)

        #Calculate fill levelfor user reference
        percent_fill =( (self.m_ox/self.V_tank) - self.rho_vap) / (self.rho_liq - self.rho_vap)
        print("\n", "ox tank % fill:", percent_fill)

        self.x_tank = ( (self.V_tank/self.m_ox) - ((self.rho_liq)**-1) )/( ((self.rho_vap)**-1) - ((self.rho_liq)**-1)) #quality
        print(self.x_tank)

        self.T_tank = CP.PropsSI('T', 'Q', self.x_tank, 'P', self.P_tank, 'N2O')  # Temperature of nitrous (K)

        self.u_tank = self.x_tank*self.u_vap + (1-self.x_tank)*self.u_liq
        self.U_tank = self.m_ox*self.u_tank

        R_UNIV = 8.314 #J/mol
        self.R = R_UNIV / CP.PropsSI('M', 'T', 300, 'P', 101325, 'N2O') #kg/mol - propsi requires dummy inputs (T, P)
        
        self.y_ox = 0

        self.rp_nos_obj = get_prop('NitrousOxide')
        self.kinematic_visc_ox = 0

        print("\n------------\nsummary of bens ox tank inputs: \nOxidizer: ", oxidizer ,"\nTimestep: ", timestep,"\nm_ox: ", m_ox ,"(kg)\nCd: ", Cd_1, "(-)\nA_inj_1: ", A_inj_1, "(m^2)\nV_tank: ", V_tank, "(m^3)\nP_tank: ", P_tank, "(Pa)\nP_cc: ", P_cc, "(Pa)\n------------\n\n\n")



    #TODO: TEST THIS FUNCTION!!!
    ### NOTE: diams input in SI!!!!!! DONT INPUT IMPERIAL FUNCTION IS CIVILIZED NOT FREE!
    def hold_time(self, T_atm, k_cond, L_wall, d_outer, d_inner):
        if(self.t > 0):
            print("WARNING HOLD TIME CALLED AT INVALID TIME T")

        v_tank = (self.V_tank/self.m_ox)

        # Get the critical temperature and pressure
        T_crit = CP.PropsSI('Tcrit', 'N2O')
        P_crit = CP.PropsSI('Pcrit', 'N2O')

        v_crit = CP.PropsSI('V', 'T', T_crit, 'P', P_crit, 'N2O')

        if(v_tank < v_crit):
            x = 1
        else:
            x = 0

        u_sat = CP.PropsSI('U', 'Q', x, 'V', v_tank, 'N2O')

        #solve Q
        Q_max = self.m_ox*(u_sat - self.u_tank)

        ####HEAT TRANSFER MODEL!!!!!!
        A_inner = 0.25*np.pi*d_inner**2
        wall_thermal_resistance = L_wall/(k_cond*A_inner)
        
        cylinder_thermal_resistance = np.ln(d_outer/d_inner)/(2*np.pi*L_wall*k_cond)

        thermal_resistance = 1/cylinder_thermal_resistance + 2/wall_thermal_resistance

        Q_dot = (T_atm - self.T_tank)/thermal_resistance

        return Q_max/Q_dot
        
          

    def inst(self, P_cc):
        self.P_cc = P_cc

        #setup iteration error tolerance
        self.V_tank_err = self.all_error*self.V_tank
        self.u_tank_err = self.all_error*self.u_tank

        ### start
        if self.x_tank < 1:
            while np.abs(Verror(self.T_tank, self.U_tank, self.m_ox, self.V_tank ) ) > self.V_tank_err:
                self.T_tank = secant((lambda T: Verror(T, self.U_tank, self.m_ox, self.V_tank)), self.T_tank)

            #use temperature to calculate thermo properties of tank
            self.P_tank = CP.PropsSI('P', 'Q', self.x_tank, 'T', self.T_tank, 'N2O')
            h_liq = CP.PropsSI('H', 'Q', 0, 'T', self.T_tank, 'N2O')
            h_vap = CP.PropsSI('H', 'Q', 1, 'T', self.T_tank, 'N2O')
            self.rho_liq = CP.PropsSI('D', 'Q', 0, 'T', self.T_tank, 'N2O')
            self.rho_vap = CP.PropsSI('D', 'Q', 1, 'T', self.T_tank, 'N2O')
            self.u_liq = CP.PropsSI('U', 'Q', 0, 'T', self.T_tank, 'N2O')
            self.u_vap = CP.PropsSI('U', 'Q', 1, 'T', self.T_tank, 'N2O')

            self.x_tank = (self.U_tank/self.m_ox - self.u_liq)/(self.u_vap - self.u_liq)
            self.u_tank = self.x_tank*self.u_vap + (1 - self.x_tank)*self.u_liq
            h_tank_exit = self.x_tank*h_vap + (1 - self.x_tank)*h_liq
            self.rho_tank = self.x_tank*self.rho_vap + (1-self.x_tank)*self.rho_liq

            #update current time
            self.t = self.t + self.timestep
            #assume only liquid draining from tank #NOTE: challenge this?
            self.rho_exit = self.rho_tank#self.rho_liq

            print("liquid phase", self.P_tank,self.x_tank, self.T_tank, self.m_ox)
            
            if(self.inj_model == 1):
                h_tank_exit = h_liq
                #SPI Model --> single phase so using liquid enthalpy
            
            

            

        else:
            print("VAPOR PHASE")
            #solve variables
            self.rho_tank = self.m_ox/self.V_tank
            self.u_tank = self.U_tank/self.m_ox

            while np.abs(uerror(self.T_tank, self.rho_tank, self.u_tank) ) > self.u_tank_err:
                self.T_tank = secant((lambda T: uerror(T, self.rho_tank, self.u_tank)), self.T_tank)

            self.P_tank = thermo_span_wagner(self.rho_tank, self.T_tank, 'p')

            h_tank_exit = thermo_span_wagner(self.rho_tank, self.T_tank, 'u')

            #update current time
            self.t = self.t + self.timestep
            self.rho_exit = self.rho_tank
            #print('here')
            h_tank_exit = h_tank_exit + 7.3397e+05 #Convert from Span-Wagner enthalpy convention to NIST


        #is this used?
        #Cp = CP.PropsSI('Cpmass', 'H', h_tank_exit, 'P', self.P_cc, 'N2O')
        #Cv = CP.PropsSI('Cvmass', 'H', h_tank_exit, 'P', self.P_cc, 'N2O')

        #self.y_ox = Cp/Cv

        ### Use Chosen Injector Model:
        #BUG: exceeding critical flow with S1 and exit density is too low which is ruining the thrust curve

        ### SPI MODEL ###
        if(self.inj_model == 1):
            #correct density for this model calculated from above based on nitrous state
            self.m_dot_ox = spi_model(self.Cd_1, self.A_inj_1, self.P_tank, self.P_cc, self.rho_exit)

            #print(self.rho_exit, self.rho_liq, self.x_tank, h_tank_exit, self.P_tank, self.P_cc)


        ### HEM MODEL ###
        elif(self.inj_model == 2):
            self.m_dot_ox = hem_model(self.Cd_1, self.A_inj_1, self.P_tank, self.P_cc, h_tank_exit)

        ### DYER MODEL ###

        #TODO: USE DYER MODEL FUNCTION ABOVE!
        elif(self.inj_model == 3):
            P_sat = CP.PropsSI('P', 'Q', 0, 'T', self.T_tank, 'N2O')
            self.m_dot_ox = dyer_model(self.Cd_1, self.A_inj_1, self.P_tank, P_sat, self.P_cc, self.rho_exit, h_tank_exit)

            #mu = CP.PropsSI('V', 'T', self.T_tank, 'P', self.P_cc, 'N2O')  # dynamic viscosity in Pa·s
            #self.kinematic_visc_ox = mu / self.rho_exit
            

        ### My implementation of Emerson and Mohammad's Modified Omega Model ###
        #https://web.stanford.edu/~cantwell/AA284A_Course_Material/AA284A_Resources/Nino%20and%20Razavi,%20Design%20of%20Two-Phase%20Injectors%20Using%20Analytical%20and%20Numerical%20Methods%20with%20Application%20to%20Hybrid%20Rockets%202019-4154.pdf
        elif(self.inj_model == 4):
            self.m_dot_ox = modified_emerson_and_mohammad_model_inst(self.P_tank, self.P_cc, self.T_tank, self.x_tank, self.A_inj_1, self.rho_exit, h_tank_exit)
        
        mu = self.rp_nos_obj.ViscAtTdegR(1.8*self.T_tank) *0.1 # convert input T from K to R and convert returned Poise to Pa s
        self.kinematic_visc_ox = mu / self.rho_exit

        #TODO: add feed system term ^

        #Ben does this to eliminate numerical instability
        if self.t == self.timestep:
            self.m_dot_ox = 0.5 * self.m_dot_ox
        else:
            self.m_dot_ox = 0.5 * self.m_dot_ox + 0.5 * self.m_dot_ox_prev

        #move forward in time with differential eqns
        print(self.x_tank, self.t)
        self.m_ox = self.m_ox - self.m_dot_ox*self.timestep
        self.m_dot_ox_prev = self.m_dot_ox
        self.U_tank = self.U_tank -self.m_dot_ox*h_tank_exit*self.timestep

        #print(self.m_dot_ox," ",self.m_ox," ",self.t)

        #print("exit rho",self.rho_exit)

        #print(self.t)