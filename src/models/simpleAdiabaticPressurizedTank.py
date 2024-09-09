###NOTE: THIS IS BAD! should have neglected mixing of helium and GOX!!!

from rocketcea.cea_obj_w_units import CEA_Obj #how to access this as a non us citizen?
import CoolProp.CoolProp as CP #I love coolprop! ~ units: http://www.coolprop.org/v4/apidoc/CoolProp.html
import matplotlib.pyplot as plt
import numpy as np

#TODO: DOUBLE CHECK PropsSI units

#TODO: FINISH!!!!!!
  
class simpleAdiabaticPressurizedTank():
    def __init__(self, pres_name, m_pres, fuel_name, m_fuel, P_tank, id_PROPTANK, TIMESTEP):
        self.pres = pres_name
        self.m_fuel = m_fuel
        self.fuel = fuel_name
        self.m_pres = m_pres
        self.P_tank = P_tank
        self.A_proptank = 0.25*np.pi*id_PROPTANK**2 #NOTE: passing in metric tank diam!!!!!

        self.TIMESTEP = TIMESTEP

        #TODO: ADD TO INIT
        self.cd = 0.4
        self.A_exit = 0.25*np.pi*(0.0254*0.5)**2
        self.T_amb = 275 #K
        self.V_tank = 0.01

        #setup!!!
        #NOTE: assuming kerosene at ambient temp!
        self.rho_prop = CP.PropsSI('D', 'P', self.P_tank, 'T', self.T_amb, fuel_name)

        V_prop = self.m_fuel/self.rho_prop
        self.T_prop = CP.PropsSI('T', 'P', self.P_tank, 'D', self.rho_prop, fuel_name)

        v_pres = (self.V_tank - V_prop)/self.m_pres

        self.T_pres = CP.PropsSI('T', 'P', self.P_tank, 'v', v_pres, pres_name)
        self.cv_pres = CP.PropsSI('CVMASS', 'P', self.P_tank, 'v', v_pres, pres_name)
        self.R_pres = 1


    def inst(self, P_downstream):

        #solve mass flow out of tank
        m_dot = self.cd * self.A_exit * np.sqrt( 2* (self.P_tank-P_downstream) / self.rho_prop )
        
        #update mass using consv of mass
        self.m_fuel -= m_dot * TIMESTEP

        #update volume, assuming incompressible kerosene
        V_fuel = self.m_fuel/self.rho_prop
        v_pres_prev = v_pres
        v_pres = (self.V_tank - V_fuel)/self.m_pres

        #solve heat transfer for conservation of energy calculation
        hc = 1
        T_helium = 1
        Q_dot = hc*self.A_proptank * (self.T_pres - self.T_prop)
        Q = Q_dot * TIMESTEP

        ###use first law to solve pressure #TODO: FILL IN HEAT TRANSFER PLACEHOLDER CONST USED FOR NOW
        self.P_tank = (self.T_pres*self.m_pres*self.cv_pres - self.P_tank*v_pres_prev +Q) / ((self.m_pres*self.cv_pres*v_pres)/(self.R_pres)+ v_pres)
        #assuming pressurant behaves as an ideal gas w const specific heats, update T w ideal gas law
        self.T_pres = self.P_tank*v_pres/self.R_pres



        #NOTE: need to solve new self.T_prop! 
        

