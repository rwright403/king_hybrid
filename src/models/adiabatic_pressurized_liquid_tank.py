###NOTE: THIS IS BAD! should have neglected mixing of helium and GOX!!!

from rocketcea.cea_obj_w_units import CEA_Obj #how to access this as a non us citizen?
import CoolProp.CoolProp as CP #I love coolprop! ~ units: http://www.coolprop.org/v4/apidoc/CoolProp.html
import matplotlib.pyplot as plt
import numpy as np


def secant(func, x1): #THANK YOU BEN KLAMMER!
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

#NOTE: not sure if this will work and if equations are too coupled, we will see
def verror(T_pres,m_pres,cv_pres,P_tank_prev,v_pres_prev,Q_transfer,v_pres,R_pres): #TODO: PASS IN ALL VARIABLES

    P_tank = (T_pres*m_pres*cv_pres + P_tank_prev*v_pres_prev + Q_transfer) / ((m_pres*cv_pres*v_pres)/(R_pres) + v_pres)
    #print("Ptank is still wrong:", P_tank, T_pres, m_pres, cv_pres, v_pres_prev, v_pres, R_pres)
        
    #assuming pressurant behaves as an ideal gas w const specific heats, update T w ideal gas law
    v_pres_new = R_pres*T_pres/P_tank

    v_tank_diff = v_pres_new - v_pres
    return v_tank_diff


#TODO: DOUBLE CHECK PropsSI units

#TODO: FINISH!!!!!!
  
class model():
    def __init__(self, pres_name, m_pres, fuel_name, m_fuel, P_tank, id_PROPTANK, V_tank_2, C_inj_2, T_amb, TIMESTEP):
        self.pres = pres_name
        self.m_fuel = m_fuel
        self.fuel = fuel_name
        self.m_pres = m_pres
        self.P_tank = P_tank
        self.id_proptank = id_PROPTANK
        self.A_proptank = 0.25*np.pi*id_PROPTANK**2 #NOTE: passing in metric tank diam!!!!!

        self.TIMESTEP = TIMESTEP

        #TODO: ADD TO INIT
        self.V_tank = V_tank_2
        self.T_amb = T_amb

        self.C_inj = C_inj_2

        #setup!!!
        #NOTE: assuming fuel density at ambient temp!
        self.rho_prop = CP.PropsSI('D', 'P', self.P_tank, 'T', self.T_amb, fuel_name)

        V_prop = self.m_fuel/self.rho_prop

        self.T_prop = CP.PropsSI('T', 'P', self.P_tank, 'D', self.rho_prop, fuel_name)

        self.v_pres = (self.V_tank - V_prop)/self.m_pres

        self.T_pres = CP.PropsSI('T', 'P', self.P_tank, 'D', (1/self.v_pres), pres_name)

        self.cv_pres = CP.PropsSI('CVMASS', 'P', self.P_tank, 'D', (1/self.v_pres), pres_name)
        self.R_pres = 8.314 / CP.PropsSI('MOLARMASS', pres_name) #8.314 J/(kg K)

        self.m_dot_fuel = 0

        ###heat transfer coefficients:
        self.C = 0.27 #NOTE: DOUBLE CHECK THIS APPLIES TO SPECIFIC FLUIDS OR IF ITS GENERAL
        self.n = 0.25
        self.K_H = 1 #this is a heat transfer corrective factor that is set to 1 from paper, might need to adjust later??
        
        self.v_tank_err = 0.000005 #NOTE: DOUBLE CHECK MAGNITUDE IF SECANT METHOD WITH VERROR FAILS WITH A SECANT ERROR

        #pres_name, m_pres, fuel_name, m_fuel, P_tank, id_PROPTANK, TIMESTEP):
        print("\n------------\nsummary of adiabatic_fuel_tank inputs: \nPressurant: ", pres_name ,"\nm_pres: ", m_pres ,"(kg)\nFuel: ", fuel_name,"\nm_fuel: ", m_fuel ,"(kg)\nP_tank: ", P_tank, "(Pa)\nid_PROPTANK: ", id_PROPTANK, "(m)\nTimestep: ", TIMESTEP,"\n------------\n\n\n")


    def inst(self, P_downstream):

        #solve mass flow out of tank
        
        self.m_dot_fuel = self.C_inj * np.sqrt( 2* self.rho_prop* (self.P_tank-P_downstream) )

        #print(self.m_dot_fuel,self.C_inj,self.P_tank,P_downstream,self.rho_prop )
        
        #update mass using consv of mass
        self.m_fuel -= self.m_dot_fuel * self.TIMESTEP

        #update volume, assuming incompressible liquid fuel
        V_fuel = self.m_fuel/self.rho_prop

        v_pres_prev = self.v_pres #TODO: DELETE IF NOT USED


        self.v_pres = (self.V_tank - V_fuel)/self.m_pres

        #solve heat transfer for conservation of energy calculation
        #https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-tankpressurization-jpp2001.pdf?emrc=66201987b6c8c
        #TODO: set this up:

        ###call coolprop to get fluid conductivity --> heat transfering from fuel to pressurant so... this should be conductivity of pressurant
        k_f = CP.PropsSI('L', 'T', self.T_pres, 'P', self.P_tank, self.pres)
        Visc_pres = CP.PropsSI('V', 'T', self.T_pres, 'P', self.P_tank, self.pres)
        cp_pres = CP.PropsSI('CPMASS', 'T', self.T_pres, 'P', self.P_tank, self.pres)
        self.cv_pres = CP.PropsSI('CVMASS', 'P', self.P_tank, 'D', (1/self.v_pres), self.pres)

        ###solve beta of pressurant:
        #solve coeff of thermal expansion
        delta_T = 1e-3  # Small temperature change for finite difference
        # Get specific volume (V) in m^3/kg
        v1 = CP.PropsSI('V', 'T', self.T_prop, 'P', self.P_tank, self.fuel)
        v2 = CP.PropsSI('V', 'T', self.T_prop - delta_T, 'P', self.P_tank, self.fuel) ###NOT SURE IF IT IS MINUS
        # Calculate the partial derivative of V with respect to T at constant P
        dVdT_P = (v2 - v1) / delta_T
        # Calculate the coefficient of thermal expansion (beta)
        self.beta = dVdT_P / v1

        #NOTE: change propellant_liq_vec to fuel and ullage to gas?
        Gr = (self.id_proptank**3)*((1/self.v_pres)**2)* 9.81 * self.beta * np.abs(self.T_prop- self.T_pres) / (Visc_pres**2) #Grashof number
        Pr = cp_pres * Visc_pres / k_f #Prandtl number

        hc = self.K_H * self.C * (k_f/self.id_proptank) * (Gr/Pr)**self.n
        Q_dot = hc*self.A_proptank * (self.T_prop - self.T_pres)
        Q_transfer= Q_dot * self.TIMESTEP


        """
        ###use first law to solve pressure #TODO: FILL IN HEAT TRANSFER PLACEHOLDER CONST USED FOR NOW
        #should be -P_tank_prev, wait no this should be ok
        #NOTE: problem is that second term is about double the first term!
        #print(v_pres_prev, self.T_pres, self.cv_pres)
        #print(self.P_tank, self.T_pres*self.m_pres*self.cv_pres, - self.P_tank*v_pres_prev, Q_transfer, ((self.m_pres*self.cv_pres*self.v_pres)/(self.R_pres)+ self.v_pres), )
        #print(self.v_pres)
        self.P_tank = (self.T_pres*self.m_pres*self.cv_pres + self.P_tank*v_pres_prev + Q_transfer) / ((self.m_pres*self.cv_pres*self.v_pres)/(self.R_pres) + self.v_pres)
        
        #v_pres_new = self.R_pres*self.T_pres/self.P_tank
        #print("testing: ",v_pres_new, self.v_pres)


        #assuming pressurant behaves as an ideal gas w const specific heats, update T w ideal gas law
        self.T_pres = self.P_tank*self.v_pres/self.R_pres

        """
        while np.abs(verror(self.T_pres,self.m_pres,self.cv_pres,self.P_tank,v_pres_prev,Q_transfer,self.v_pres,self.R_pres) ) > self.v_tank_err:
            print("looping", self.T_pres)
            self.T_pres = secant((lambda T: verror(T, self.m_pres,self.cv_pres,self.P_tank,v_pres_prev,Q_transfer,self.v_pres,self.R_pres)), self.T_pres)
        
        #now use new T_pres to solve P_tank
        self.P_tank = self.R_pres*self.T_pres/self.v_pres
        
        #print(self.v_pres, self.P_tank, self.T_pres)

        #print("Ptank is still wrong:", self.P_tank, self.T_pres, self.m_pres, self.cv_pres, v_pres_prev, self.v_pres, self.R_pres)
        #print(self.P_tank, self.m_fuel, self.v_pres, (self.V_tank - V_fuel), "m^3")


        

