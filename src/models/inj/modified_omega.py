import numpy as np
from scipy.optimize import root_scalar
import CoolProp.CoolProp as CP
from src.models.inj._base import BaseInjector
from src.models.inj.dyer import dyer_model


def TWOPHASEerror(eta_crit, omega):
    function_diff = (eta_crit**2) + ((omega**2)-2*omega)*((1-eta_crit)**2) + 2*(omega**2)*np.log(eta_crit) + 2*(omega**2)*(1-eta_crit)  
    return function_diff

def LOWSUBCOOLEDerror(eta_crit, eta_sat, omega_sat):
    function_diff = (((omega_sat+(1/omega_sat)-2)/(2*eta_sat))*(eta_crit**2)) - (2*(omega_sat-1)*eta_crit) + (omega_sat*eta_sat*np.log(eta_crit/eta_sat)) + ((3/2)*omega_sat*eta_sat) - 1
    return function_diff


class modified_omega_model(BaseInjector):
    """
    My implementation of Emerson and Mohammad's Modified Omega Model
    https://web.stanford.edu/~cantwell/AA284A_Course_Material/AA284A_Resources/Nino%20and%20Razavi,%20Design%20of%20Two-Phase%20Injectors%20Using%20Analytical%20and%20Numerical%20Methods%20with%20Application%20to%20Hybrid%20Rockets%202019-4154.pdf
    """

    def __init__(self, Cd: float, A_inj: float):
        super().__init__(Cd, A_inj)
        self.A_inj = A_inj
        self.dyer = dyer_model(Cd, A_inj)
        self.fluid = "N2O"



    def m_dot(self, state: dict) -> float:
        P_2, P_1, T_1, x_1, rho_1, h_1 = state["P_2"], state["P_1"], state["T_1"], state["x_1"], state["rho_1"], state["h_1"]

        all_err = 0.005 #TODO: maybe UPDATE TO USE WHAT THE REST OF THE FUNCTION USES?


        #1) --> check inlet fluid state (sat liq. or subcooled)
        phase = CP.PropsSI('Phase', 'P', P_1, 'D', rho_1, 'N2O')




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
        
        # eta_sat is the initial guess for critical pressure ratio

        sol = root_scalar(
            lambda eta: TWOPHASEerror(eta, omega_sat),
            method="secant",
            x0=eta_sat,
            x1=eta_sat * 0.99,
            xtol=all_err,
            maxiter=100
        )
        if not sol.converged:
            raise RuntimeError("TWOPHASEerror solver did not converge")
        eta_crit_sat = sol.root


        #while np.abs(TWOPHASEerror(eta_crit_sat, omega_sat) ) > all_err:
        #    eta_crit_sat = secant((lambda T: TWOPHASEerror(T, omega_sat)), eta_crit_sat)
        

        #print("eta sat: ", eta_sat, eta_crit_sat)
        P_crit_sat = eta_crit_sat * P_1
        #print("\n", eta_crit_sat * P_1, eta_crit_sat * P_sat, "\n")

        if P_crit_sat >= P_2: #Saturated Inlet Choked
            G_sat = eta_crit_sat * np.sqrt(P_1*rho_1/omega_sat)

            print("sat inlet choked")

        else: #Saturated Inlet Not Choked
            print("sat inlet not choked")
            pratio = P_2/P_sat
            G_sat = np.sqrt(P_1*rho_1) *np.sqrt( -2*(omega_sat * np.log(pratio) + (omega_sat -1)*(1 - pratio)) ) / (omega_sat*((1/(pratio)) - 1) + 1)



            #print("sat inlet not choked", G_sat,  (-2*(omega_sat * np.log(eta_sat))) , ((omega_sat -1)*(1 - eta_sat)), omega_sat, eta_sat)

        if phase == 6: #saturated liquid vapor
            k_cavitation_const = (P_1 - P_sat) / (P_1 - P_2)
            if k_cavitation_const == 0:
                k_cavitation_const = 1 #cavitation constant only valid for P_sat > P_2, otherwise K = 1
            Cd_twophase = 0.386 + 0.361*np.sqrt(k_cavitation_const)
            m_dot =  ( Cd_twophase * self.A_inj) * G_sat  
            print("sat liq vap")

            
            #print("line 232", m_dot, rho_1, P_1, eta_crit_sat, omega_sat, x_1)
            #checked A_inj_1 is correct, fixed issue with cavitation constant k

        
        elif phase == 0: #subcooled fluid 

            #setup subcooled fluid inlet 
            
            rho_1_l = CP.PropsSI('D', 'Q', 0, 'P', P_1, 'N2O')
            eta_transition =  2*omega_sat / (1 + 2*omega_sat)

            ### High subcooled
            if P_sat <= (eta_transition * P_1):
                print("HIGH SUBCOOLED")

                P_sat = 0.9*CP.PropsSI('P', 'T', T_1, 'Q', 0, 'N2O')
                #print("correction factor of 0.9")
                omega_sat = (c_1_l*T_1*P_sat/v_1_l)*( (v_1_lg/h_1_lg)**2) #TODO: put corr factor on omega instead?
                eta_transition =  2*omega_sat / (1 + 2*omega_sat)

                eta_crit = (P_sat / (P_1)) 
                P_crit = P_sat

                if P_2 < P_crit: #choked flow
                    Cd_high_supercharge = 0.73
                    m_dot = (Cd_high_supercharge*self.A_inj) * np.sqrt( 2*(1-eta_crit)*P_1*rho_1)

                else: 
                    m_dot = 1.105*dyer_model( Cd_hem_spi_dyer, self.A_inj, P_1, P_sat, P_2, rho_1, h_1 )




            ### Low subcooled
            else:
                print("low subcooled")
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

                #while np.abs(LOWSUBCOOLEDerror(eta_crit_low, eta_sat, omega_sat) ) > all_err:
                #    eta_crit_low = secant((lambda T: LOWSUBCOOLEDerror(T, eta_sat, omega_sat)), eta_crit_low)
                
                sol = root_scalar(
                    lambda eta: LOWSUBCOOLEDerror(eta, eta_sat, omega_sat),
                    method="secant",
                    x0=eta,
                    x1=eta * 0.99,
                    xtol=all_err,
                    maxiter=100
                )
                if not sol.converged:
                    raise RuntimeError("LOWSUBCOOLEDerror solver did not converge")
                eta_crit_low = sol.root


                P_crit_low = eta_crit_low * P_1

                P_crit = eta_sat*P_crit_sat + (1-eta_sat)*P_crit_low #NOTE: using modified omega model to predict choking pressure and mass flow rate

                
                #check for choking
                if P_2 <= P_crit: #choked flow
                    G_low = np.sqrt(rho_1_l * P_1) * np.sqrt( 2*(1-eta_sat) + 2*(omega_sat*eta_sat*np.log(eta_sat/eta_crit_low) - (omega_sat-1)*(eta_sat-eta_crit_low))) / (omega_sat*((eta_sat/eta_crit_low) - 1) + 1)
                    m_dot= self.A_inj *( (P_sat/P_1)*G_sat + (1-(P_sat/P_1))*G_low )

                else: #not choked use dyer model --> smoothing factor based on choked mass flow rate!

                    #solve choked: 
                    G_sat_choked = eta_crit_sat * np.sqrt(P_1*rho_1/omega_sat)
                    #implicitly solve eta_crit_low
                    eta = P_2 / P_1
                    #eta_crit_low = eta #initial guess for critical pressure ratio

                    sol = root_scalar(
                        lambda eta: LOWSUBCOOLEDerror(eta, eta_sat, omega_sat),
                        method="secant",
                        x0=eta,
                        x1=eta * 0.99,
                        xtol=all_err,
                        maxiter=100
                    )
                    if not sol.converged:
                        raise RuntimeError("LOWSUBCOOLEDerror solver did not converge")
                    eta_crit_low = sol.root

                    #while np.abs(LOWSUBCOOLEDerror(eta_crit_low, eta_sat, omega_sat) ) > all_err:
                    #    eta_crit_low = secant((lambda T: LOWSUBCOOLEDerror(T, eta_sat, omega_sat)), eta_crit_low)
                    
                    G_low_choked = np.sqrt(rho_1_l * P_1) * np.sqrt( 2*(1-eta_sat) + 2*(omega_sat*eta_sat*np.log(eta_sat/eta_crit_low) - (omega_sat-1)*(eta_sat-eta_crit_low))) / (omega_sat*((eta_sat/eta_crit_low) - 1) + 1)
                    m_dot_choked = self.A_inj *( (P_sat/P_1)*G_sat_choked + (1-(P_sat/P_1))*G_low_choked )
                    
                    smoothing_factor = (6.083086 * m_dot_choked + 0.8304554896) #I made this up to get the transition to be smooth
                    #print("smoothing factor: ", smoothing_factor)
                    m_dot = (smoothing_factor)*dyer_model( Cd_hem_spi_dyer, self.A_inj, P_1, P_sat, P_2, rho_1, h_1 )

        return m_dot
