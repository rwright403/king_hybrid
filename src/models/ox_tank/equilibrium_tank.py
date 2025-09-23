import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import root_scalar
from src.models.inj._base import build_state
from src.models._thermo.n2o_thermo_span_wagner_class import SpanWagnerEOS_SingleState
from src.models.ox_tank._base import BaseTank


def Verror(T_tank, U_tank, m_ox, V_tank):

    rho_liq = CP.PropsSI('D', 'Q', 0, 'T', T_tank, 'N2O')  # Density of liquid nitrous (kg/m^3)
    rho_vap = CP.PropsSI('D', 'Q', 1, 'T', T_tank, 'N2O')  # Density of nitrous gas (kg/m^3)

    u_liq = CP.PropsSI('U', 'Q', 0, 'T', T_tank, 'N2O')  # Internal Energy of liquid nitrous (J/kg)
    u_vap = CP.PropsSI('U', 'Q', 1, 'T', T_tank, 'N2O')  # Internal Energy of nitrous gas (J/kg)

    x_tank = (U_tank/m_ox - u_liq)/(u_vap - u_liq)
    V_tank_diff = m_ox * ( ((1-x_tank) /rho_liq) + (x_tank/rho_vap) ) - V_tank
    return V_tank_diff

def uerror(T, rho_tank, u_tank):
    u2 = SpanWagnerEOS_SingleState(rho_tank, T).u
    u_diff = u2 - u_tank
    return u_diff


"""
Adapted From Ben Klammer's Mech 498 work
"""
class equilibrium_tank_model(BaseTank):
    def __init__(self, timestep, m_ox, V_tank, P_tank, P_atm, all_error, injector):
            super().__init__(injector, timestep)
            self.oxidizer = "N2O"
            self.m_ox = m_ox
            self.V_tank = V_tank
            self.P_tank = P_tank
            self.all_error = all_error


            self.timestep = timestep
            self.P_cc = P_atm
            self.all_error = all_error




            #setup initial tank condition
            self.t = 0

            self.m_dot_ox = 0
            self.m_dot_ox_prev = self.m_dot_ox

            #setup thermo data
            self.rho_liq = CP.PropsSI('D', 'Q', 0, 'P', self.P_tank, 'N2O')  # Density of liquid nitrous (kg/m^3)
            self.rho_vap = CP.PropsSI('D', 'Q', 1, 'P', self.P_tank, 'N2O')  # Density of nitrous gas (kg/m^3)

            self.u_liq = CP.PropsSI('U', 'Q', 0, 'P', self.P_tank, 'N2O')  # Internal Energy of liquid nitrous (kJ/kg)
            self.u_vap = CP.PropsSI('U', 'Q', 1, 'P', self.P_tank, 'N2O')  # Internal Energy of nitrous gas (kJ/kg)

            self.x_tank = ( (self.V_tank/self.m_ox) - ((self.rho_liq)**-1) )/( ((self.rho_vap)**-1) - ((self.rho_liq)**-1)) #quality

            self.T_tank = CP.PropsSI('T', 'Q', self.x_tank, 'P', self.P_tank, 'N2O')  # Temperature of nitrous (K)

            self.u_tank = self.x_tank*self.u_vap + (1-self.x_tank)*self.u_liq
            self.U_tank = self.m_ox*self.u_tank

            self.y_ox = 0

            #setup iteration error tolerance
            self.V_tank_err = self.all_error*self.V_tank
            self.u_tank_err = self.all_error*self.u_tank


            #print("\n------------\nsummary of bens ox tank inputs: \nOxidizer: ", oxidizer ,"\nTimestep: ", timestep,"\nm_ox: ", m_ox ,"(kg)\nCd: ", Cd_1, "(-)\nA_inj_1: ", A_inj_1, "(m^2)\nV_tank: ", V_tank, "(m^3)\nP_tank: ", P_tank, "(Pa)\nP_cc: ", P_cc, "(Pa)\n------------\n\n\n")     
            #Calculate fill levelfor user reference
            #percent_fill =( (self.m_ox/self.V_tank) - self.rho_vap) / (self.rho_liq - self.rho_vap)
            #print("\n", "ox tank % fill:", percent_fill)

                        # Build state dict for injector
            self.state = build_state()

    




    def inst(self, P_cc: float):


        ### start

        #print("x_tank! ", self.x_tank)
        if self.x_tank < 1:

            sol = root_scalar(
                lambda T: Verror(T, self.U_tank, self.m_ox, self.V_tank),
                method="secant",
                x0=self.T_tank,
                x1=self.T_tank * 0.99,
                xtol=self.all_error,
                maxiter=1000
            )
            if not sol.converged:
                raise RuntimeError("root_scalar failed to converge in equilibrium tank (liquid phase)")
            self.T_tank = sol.root


            self.P_tank = CP.PropsSI('P', 'Q', self.x_tank, 'T', self.T_tank, 'N2O')

            h_liq = CP.PropsSI('H', 'Q', 0, 'T', self.T_tank, 'N2O')
            h_vap = CP.PropsSI('H', 'Q', 1, 'T', self.T_tank, 'N2O')

            self.rho_liq = CP.PropsSI('D', 'Q', 0, 'T', self.T_tank, 'N2O')
            self.rho_vap = CP.PropsSI('D', 'Q', 1, 'T', self.T_tank, 'N2O')

            self.u_liq = CP.PropsSI('U', 'Q', 0, 'T', self.T_tank, 'N2O')
            self.u_vap = CP.PropsSI('U', 'Q', 1, 'T', self.T_tank, 'N2O')

            self.x_tank = (self.U_tank/self.m_ox - self.u_liq)/(self.u_vap - self.u_liq)

            self.u_tank = (self.x_tank*self.u_vap + (1 - self.x_tank)*self.u_liq )

            h_tank_exit = self.x_tank*h_vap + (1 - self.x_tank)*h_liq

            self.rho_tank = self.x_tank*self.rho_vap + (1-self.x_tank)*self.rho_liq


            #assume only liquid draining from tank
            self.rho_exit = self.rho_liq

            #print("liquid phase", self.P_tank,self.P_cc,self.x_tank,self.rho_tank,self.rho_liq, self.t)
            self.state["P_sat"] = self.P_tank
                
            

            

        else:
            print("VAPOR PHASE")
            self.rho_tank = self.m_ox/self.V_tank
            self.u_tank = self.U_tank/self.m_ox

            sol = root_scalar(
                lambda T: uerror(T, self.rho_tank, self.u_tank),
                method="secant",
                x0=self.T_tank,
                x1=self.T_tank * 0.99,
                xtol=self.all_error,
                maxiter=1000
            )
            if not sol.converged:
                raise RuntimeError("root_scalar failed to converge in equilibrium tank (vapor phase)")
            self.T_tank = sol.root

            self.P_tank = SpanWagnerEOS_SingleState(self.rho_tank, self.T_tank).P

            h_tank_exit = SpanWagnerEOS_SingleState(self.rho_tank, self.T_tank).u

            self.rho_exit = self.rho_tank
            h_tank_exit = h_tank_exit



        self.state["P_1"] = self.P_tank
        self.state["P_2"] = P_cc
        self.state["rho_1"] = self.rho_exit
        self.state["h_1"] = h_tank_exit
        self.state["T_1"] = self.T_tank
        self.state["x_1"] = self.x_tank

        self.m_dot_ox = self.injector.m_dot(self.state)

        self.t = self.t + self.timestep #update current time

        #Ben does this to eliminate numerical instability
        
        if self.t == self.timestep:
            self.m_dot_ox = 0.5 * self.m_dot_ox
        else:
            self.m_dot_ox = 0.5 * self.m_dot_ox + 0.5 * self.m_dot_ox_prev
        

        self.m_ox = self.m_ox - self.m_dot_ox*self.timestep
        self.m_dot_ox_prev = self.m_dot_ox
        self.U_tank = self.U_tank -self.m_dot_ox*h_tank_exit*self.timestep

        return {"P_ox_tank": self.P_tank, "m_dot_ox_tank": self.m_dot_ox, "m_ox_tank": self.m_ox, "T_ox_tank": self.T_tank}


