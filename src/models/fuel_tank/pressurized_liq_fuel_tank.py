import numpy as np
from src.utils.numerical_methods import rk4_step

import CoolProp.CoolProp as CP

from src.utils.enum import FillType
from src.models.inj._base import build_state
from src.models.inj.spi import spi_model
from src.models._thermo.convection_heat_transfer import *
from src.models._thermo.conduction_heat_transfer import *
from src.models.fuel_tank._base import BaseTank

#TODO: MAKE AN INPUT
CW = 896 #J/(kg K) ~ this is the specific heat capacity of the wall material, for Al 6061 from MATWEB: https://www.matweb.com/search/datasheet.aspx?MatGUID=b8d536e0b9b54bd7b69e4124d8f1d20a&ckck=1
    
#TODO: NEED TO FIND FUEL HEAT TRANSFER COEFFS!!!!

def adiabatic_fill(m_total, V_pres_total, U_pres_total, pres_str):
    """
    Assume adiabatic fill, total energy between valve openning and fuel tank being pressurized is conserved.
    This holds best assuming a fast process and we launch quickly after pressurizing our fuel tank
    """

    u_adiabatic_fill = U_pres_total/m_total
    rho_adiabatic_fill = m_total/V_pres_total

    P_tank = CP.PropsSI('P', 'U', u_adiabatic_fill, 'D', rho_adiabatic_fill, pres_str)

    return rho_adiabatic_fill, P_tank



def isothermal_fill(m_total, V_pres_total, T_atm, pres_str):
    """
    Assume isothermal fill, temperature of pressurized gas reaches eq with env before launch.
    This holds best assuming we pressurize and sit on the pad for awhile before launch.
    """

    rho_isothermal_fill = m_total/V_pres_total

    P_tank = CP.PropsSI('P', 'T', T_atm, 'D', rho_isothermal_fill, pres_str)

    return rho_isothermal_fill, P_tank



class pressurized_liq_fuel_tank(BaseTank):
    def __init__(self, timestep, m_fuel, m_pres, P_pres_tank, P_atm, T_atm, rho_atm, V_fuel_tank, V_pres_tank, diam_out, diam_in, rho_wall, k_w, Cd, A_inj, fuel_str, pres_str, fill_type:FillType ):
        super().__init__(injector =spi_model(Cd, A_inj), timestep=timestep)
        #self.spi_model(Cd, A_inj) #since assume only liquid draining from tank, we are forcing this model to be the SPI model
        self.state = build_state()

        self.timestep = timestep

        # Environment 
        self.T_atm = T_atm
        self.P_atm = P_atm
        self.rho_atm = rho_atm

        self.fuel_str = fuel_str
        self.pres_str = pres_str


        # Fuel Tank Geometry
        if diam_in >= diam_out:
            raise ValueError("Tank inner diameter >= Tank outer diameter, this will give nonphysical results so check your inputs!")

        self.P_tank = P_pres_tank
        self.V_fuel_tank = V_fuel_tank

        self.diam_out = diam_out
        self.diam_in = diam_in
        self.Tank_Inner_Area = 0.25*np.pi*(self.diam_in**2)
        self.Tank_xsec_Area = 0.25*np.pi*((self.diam_out**2)-(self.diam_in**2))
        self.height_tank = self.V_fuel_tank/self.Tank_Inner_Area #resolve tank height for new volume

        self.rho_wall = rho_wall
        self.k_w = k_w

        self.m_fuel = m_fuel
        self.rho_fuel = CP.PropsSI('D', 'P', P_pres_tank, 'T', T_atm, self.fuel_str)
        self.state["rho_1"] = self.rho_fuel # incompressible fuel assumption

        V_fuel = self.m_fuel/self.rho_fuel
        self.T_fuel = CP.PropsSI('T', 'P', self.P_tank, 'D', self.rho_fuel, self.fuel_str)

        V_pres_ullage = self.V_fuel_tank - V_fuel #Problem, this is (-)
        self.T_pres = self.T_atm
        self.rho_pres = CP.PropsSI('D', 'P', self.P_tank, 'T', self.T_pres, self.pres_str) 
        self.m_pres = V_pres_ullage*self.rho_pres



        self.T_wall_fuel = T_atm
        self.T_wall_pres = T_atm
        

        print("init fuel tank: ", self.m_pres, self.rho_pres, self.T_pres)

        """
        Setup Initial fuel_tank. 
        Start by initializing pressurant tank and opening valve 
        to pressurize fuel_tank before first simulated timestep!
        Assuming:
        - Fuel in tank starts at T_atm, P_atm
        - Air in tank at start is treated as pressurant gas,
          which is also at T_atm and P_atm
        - Pressurant tank and Ullage in fuel tank form same ctrl volume
        - Pressurant is NOT soluble in fuel, so we dont need to worry about that problem
        """
        """
        # Sol mass of fuel at T_atm, P_atm in ullage space
        P_tank = P_atm
        # Fuel assumed to be incompressible!
        self.rho_fuel = CP.PropsSI('D', 'P', P_tank, 'T', T_atm, self.fuel_str)
        self.m_fuel = m_fuel
        self.V_fuel = self.m_fuel/self.rho_fuel
        V_pres_ullage = self.V_fuel_tank - self.V_fuel 

        #print("initial conditions: ", self.rho_fuel, self.m_fuel, self.V_fuel, V_pres_ullage)

        
        rho_pres = CP.PropsSI('D', 'P', P_tank, 'T', T_atm, self.pres_str)

        m_pres_ullage = rho_pres*V_pres_ullage
        u_pres_tank = CP.PropsSI('U', 'P', P_pres_tank, 'T', T_atm, self.pres_str)
        U_pres_tank = u_pres_tank*m_pres

        u_pres_ullage = CP.PropsSI('U', 'P', self.P_atm, 'T', T_atm, self.pres_str)
        U_pres_ullage = u_pres_ullage*m_pres_ullage

        U_pres_total = U_pres_tank + U_pres_ullage
        V_pres_total = V_pres_tank + V_pres_ullage
        m_total = m_pres + m_pres_ullage

        self.rho_pres, self.P_tank = None, None

        if fill_type == FillType.ADIABATIC:
            self.rho_pres, self.P_tank = adiabatic_fill(m_total, V_pres_total, U_pres_total, self.pres_str)
        elif fill_type == FillType.ISOTHERMAL:
            self.rho_pres, self.P_tank = isothermal_fill(m_total, V_pres_total, T_atm, self.pres_str)
        else:
            raise ValueError(f"Unknown fill type: {fill_type}")

        # Now resolve system as if it was just one tank

        self.T_wall_fuel = T_atm
        self.T_wall_pres = T_atm

        self.m_pres = m_total 
        self.V_fuel_tank += V_pres_total #this is messy fix
        self.height_tank = self.V_fuel_tank/self.Tank_Inner_Area #resolve tank height for new volume


        self.T_fuel = CP.PropsSI('T', 'P', self.P_tank, 'D', self.rho_fuel, self.fuel_str)
        self.T_pres = CP.PropsSI('T', 'P', self.P_tank, 'D', self.rho_pres, self.pres_str) 
        self.state["rho_1"] = self.rho_fuel # incompressible fuel assumption
        """

    def system_of_fuel_pres_odes(self, t, y, P_cc):               

        T_fuel, T_pres, m_fuel, T_wall_fuel, T_wall_pres, _ = y  # Unpack state variables

        ### Solve thermo parameters!
        V_fuel = m_fuel/self.rho_fuel # incompressible fuel assumption
        V_pres = self.V_fuel_tank - V_fuel
        rho_pres = self.m_pres/V_pres

        P_tank = CP.PropsSI('P', 'D', rho_pres, 'T', T_pres, self.pres_str)

        #print("pres liq fuel init: ", P_tank, self.rho_fuel, rho_pres )

        # fuel
        T_fuel = CP.PropsSI('T', 'P', P_tank, 'D', self.rho_fuel, self.fuel_str)
        h_fuel = CP.PropsSI('HMASS', 'T', T_fuel, 'P', P_tank, self.fuel_str)
        u_fuel = CP.PropsSI('UMASS', 'T', T_fuel, 'P', P_tank, self.fuel_str)
        cv_fuel = CP.PropsSI('CVMASS', 'T', T_fuel, 'P', P_tank, self.fuel_str)

        # pressurant
        cv_pres = CP.PropsSI('CVMASS', 'T', T_pres, 'D', rho_pres, self.pres_str)

        ### Solve Mass Balance
        self.state["P_1"] = self.P_tank
        self.state["P_2"] = P_cc
        m_dot_fuel = -self.injector.m_dot(self.state)

        ### Propellant node heat transfer:

        #then solve the height of the pres wall
        h_pres_wall = V_pres / self.Tank_Inner_Area
        V_pres_wall = self.Tank_xsec_Area*h_pres_wall
        m_pres_wall = self.rho_wall*V_pres_wall

        h_fuel_wall = self.height_tank - h_pres_wall
        V_fuel_wall = self.Tank_xsec_Area*h_fuel_wall
        m_fuel_wall = self.rho_wall*V_fuel_wall

        #print("pres fuel tank into heat transfer: ", self.rho_fuel, rho_pres, T_wall_fuel, T_fuel, T_fuel, P_tank, h_fuel_wall, (np.pi*self.diam_in*h_fuel_wall + self.Tank_Inner_Area))
        #print("h_fuel_wall should not be negative! ", h_fuel_wall)


        # Heat transfer (4) [natural convection] from fuel wall to fuel
        Q_dot_fuel_wall_to_fuel = 0#solve_Q_dot_natural_convection_liq(self.rho_fuel, T_wall_fuel, T_fuel, T_fuel, P_tank, 0.021, 0.4, h_fuel_wall, (np.pi*self.diam_in*h_fuel_wall + self.Tank_Inner_Area), self.fuel_str ) #relative to fuel cv       
        # Heat transfer (5) [natural convection] from pres to pres wall
        Q_dot_pres_wall_to_pres = 0#solve_Q_dot_natural_convection_gas(rho_pres, T_wall_pres, T_pres, T_pres, P_tank, 0.021, 0.4, h_pres_wall, (np.pi*self.diam_in*h_pres_wall + self.Tank_Inner_Area), self.pres_str ) #relative to pres cv

        Q_dot_fuel_to_pres = solve_Q_dot_natural_convection_gas(rho_pres, T_fuel, T_pres, T_pres, P_tank, 0.021, 0.4, h_pres_wall, (np.pi*self.diam_in*h_pres_wall + self.Tank_Inner_Area), self.pres_str ) #relative to pres cv


        Q_dot_fuel = Q_dot_fuel_wall_to_fuel - Q_dot_fuel_to_pres
        Q_dot_pres = Q_dot_pres_wall_to_pres + Q_dot_fuel_to_pres


        ### Wall nodes:

        #NOTE: USE ambient properties for air, T_2 will be respective wall temperature (RK var)
        # (6) [natural convection] from atm to liq wall
        Q_dot_atm_to_fuel_wall = 0#solve_Q_dot_natural_convection_gas(self.rho_atm, self.T_atm, T_wall_fuel, self.T_atm, self.P_atm, 0.59, 0.25, self.height_tank, (np.pi*self.diam_out*h_fuel_wall + self.Tank_Inner_Area), "Air") #relative to wall_liq cv
        # (7) [natural convection] from atm to gas wall
        Q_dot_atm_to_pres_wall = 0#solve_Q_dot_natural_convection_gas(self.rho_atm, self.T_atm, T_wall_pres, self.T_atm, self.P_atm, 0.59, 0.25, self.height_tank, (np.pi*self.diam_out*h_pres_wall + self.Tank_Inner_Area), "Air") #relative to wall_gas cv
        # (8) [conduction] from liq wall to gas wall 
        Q_dot_fuel_wall_to_pres_wall = 0#solve_Q_dot_conduction( (T_wall_fuel-T_wall_pres), self.height_tank, self.k_w, self.Tank_xsec_Area) #relative to wall_liq cv


        ### solve fuel height rate of change
        V_dot_fuel = m_dot_fuel/self.rho_fuel
        height_dot = V_dot_fuel/self.Tank_Inner_Area

        m_dot_liq_wall = self.rho_wall*(0.25*np.pi*height_dot*((self.diam_out**2)-(self.diam_in**2)))  #BUG: this might be a bit unstable w runge kutta steps?
        m_dot_gas_wall = -m_dot_liq_wall

        print("V_dot_fuel: ", V_dot_fuel)

        ### solve 
        #NOTE: for T_dot_wall_liq:  IN: (6)      OUT: (4) AND (8)
        T_dot_wall_fuel = ( Q_dot_atm_to_fuel_wall - Q_dot_fuel_wall_to_fuel - Q_dot_fuel_wall_to_pres_wall + m_dot_liq_wall*CW*(T_wall_fuel - T_wall_pres) ) / (CW*m_fuel_wall)

        #NOTE: for T_dot_wall_vap:  IN: (7) and (8)      OUT: (5)
        T_dot_wall_pres = ( Q_dot_atm_to_pres_wall - Q_dot_pres_wall_to_pres + Q_dot_fuel_wall_to_pres_wall + m_dot_gas_wall*CW*(T_wall_fuel - T_wall_pres) ) / (CW*m_pres_wall)


        ### Solve Propellant Nodes Energy Balance

        U_dot_fuel = m_dot_fuel*h_fuel - P_tank*V_dot_fuel + Q_dot_fuel
        U_dot_pres = (-1)*P_tank*(-1*V_dot_fuel) + Q_dot_pres #NOTE: V_dot_pres = (-1)*V_dot_fuel

        print("U_dot_pres: ", U_dot_pres , (-1)*P_tank*(-1*V_dot_fuel) , Q_dot_pres )


        T_dot_fuel = (U_dot_fuel - u_fuel*m_dot_fuel)/(m_fuel*cv_fuel)


        rho_dot_pres = -(self.m_pres/(V_pres**2))*(-1*V_dot_fuel) #NOTE: V_dot_pres = (-1)*V_dot_fuel
        
        du_drho_const_T_pres = CP.PropsSI('d(UMASS)/d(Dmass)|T', 'T',  T_pres, 'D', rho_pres, self.pres_str)

        T_dot_pres = (1/cv_pres)*( (1/self.m_pres) * (U_dot_pres) - (du_drho_const_T_pres * rho_dot_pres) )


        print("T_dot_pres: ", T_dot_pres, (1/self.m_pres) * (U_dot_pres), - (du_drho_const_T_pres * rho_dot_pres) )

        return [T_dot_fuel, T_dot_pres, m_dot_fuel, T_dot_wall_fuel, T_dot_wall_pres, rho_dot_pres]

    """

    def system_of_inert_ullage_pres_odes(self, t, y, P_cc):

        # Unpack state
        _, T_pres, _, _, T_wall_pres, rho_pres = y


        # Convert density -> mass (for energy ODE)
        m_pres = max(rho_pres * self.V_tank, 0.0)

        # Tank pressure from ideal gas
        self.P_tank = CP.PropsSI('P', 'D', rho_pres, 'T', T_pres, self.pres_str)

        # --- orifice mass flow (reservoir -> back pressure P_cc) ---
        # guard against non-physical states
        if P_tank <= P_cc or P_tank <= 0.0 or T_pres <= 0.0 or m_pres <= 0.0:
            mdot_out = 0.0
        else:
            pr = P_cc / P_tank
            pr_crit = (2.0/(g+1.0))**(g/(g-1.0))
            coeff = CdA * P_tank * (g/(R*T_pres))**0.5



            if pr <= pr_crit:  # choked
                phi_star = (2.0/(g+1.0))**((g+1.0)/(2.0*(g-1.0)))
                mdot_out = coeff * phi_star
            else:              # unchoked
                term = (2.0/(g-1.0)) * (1.0 - pr**((g-1.0)/g))
                mdot_out = coeff * (pr**(1.0/g)) * max(term, 0.0)**0.5

        # --- ODEs (adiabatic, rigid tank) ---
        # Mass: dm/dt = - m_dot_out      ->  d(rho)/dt = (1/V) dm/dt
        dm_pres_dt   = -mdot_out
        drho_pres_dt = dm_pres_dt / self.V_tank

        # Energy in tank:  m c_v dT/dt = - m_dot_out * h_out  with h_out = c_p * T
        # Using c_p = gamma*R/(gamma-1), c_v = R/(gamma-1)  => dT/dt = -(gamma-1)*(m_dot/m)*T
        if m_pres > 0.0:
            dT_pres_dt = -(g - 1.0) * (mdot_out / m_pres) * T_pres
        else:
            dT_pres_dt = 0.0

        # --- everything else "off" per your request ---
        dT_fuel_dt       = 0.0
        dm_fuel_dt       = 0.0
        dT_wall_fuel_dt  = 0.0
        dT_wall_pres_dt  = 0.0



        return [0.0, dT_pres_dt, 0.0, 0.0, 0.0, drho_pres_dt]


    """



    def tank_ode_system(self, t, y, P_cc):
        """
        Wrapper that chooses liquid-phase or vapor-only ODEs.
        """
        _, _, m_fuel, _, _, _ = y

        if m_fuel > 0:
            return self.system_of_fuel_pres_odes(t, y, P_cc)
        else:
            # fuel fully drained from tank
            print("gas phase!")

            return [self.T_fuel, self.T_pres, self.m_fuel, self.T_wall_fuel, self.T_wall_pres, 0.0] #NOTE: rho_dot_pres = 0.0 (density not changing just holding tank const)


    def inst(self, P_cc):

        y0 = [self.T_fuel, self.T_pres, self.m_fuel, self.T_wall_fuel, self.T_wall_pres, self.rho_pres]
        y_new = rk4_step(self.tank_ode_system, 0.0, y0, self.timestep, P_cc)

        #NOTE: this only applies if liq phase: 
        if self.m_fuel > 0:
            self.T_fuel, self.T_pres, self.m_fuel, self.T_wall_fuel, self.T_wall_pres, self.rho_pres = y_new

            self.P_tank = CP.PropsSI('P', 'D', self.rho_pres, 'T', self.T_pres, self.pres_str)


        self.state["P_1"] = self.P_tank
        self.state["P_2"] = P_cc
        m_dot_fuel = self.injector.m_dot(self.state)

        return {"P_fuel_tank": self.P_tank, "m_dot_fuel_tank": m_dot_fuel, "m_fuel": self.m_fuel}
    
    #TODO: NEED TO VENT ULLAGE GAS WITHOUT RETURNING "FUEL" MASS FLOW RATE TO CC