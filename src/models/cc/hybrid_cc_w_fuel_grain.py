#https://rocketcea.readthedocs.io/en/latest/functions.html#module-rocketcea.cea_obj
import numpy as np
from scipy.optimize import root_scalar
from rocketcea.cea_obj import add_new_fuel
from rocketcea.cea_obj_w_units import CEA_Obj
from src.models.cc._base import BaseChamber
from src.utils.numerical_methods import rk4_step

R_UNIV = 8314  # J / (kmol K)



import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
from src.models.cc._base import BaseChamber

R_UNIV = 8314.0  # J/kmol-K


class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, L_star, m_fuel_i, rho_fuel, a, n, L, A_port, nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        self.C = C
        self.rho_fuel = rho_fuel
        self.m_fuel_t = m_fuel_i
        self.a = a
        self.n = n

        # geometry
        self.L = L
        self.A_port = A_port
        self.radius = np.sqrt(A_port / np.pi)

        # nozzle
        self.expratio = nozzle.A_exit / nozzle.A_throat

        # timestep
        self.timestep = timestep

        # chamber state
        self.P_cc = P_atm
        self.T_cc = 298.0
        self.R = 0.0
        self.gamma = 1.2
        self.cp = 1000.0

        # flow + regression
        self.r_dot_t = 0.0
        self.m_dot_fuel = 0.0
        self.m_dot_cc_t = 0.0
        self.OF = 0.0

        # nozzle + thrust
        self.v_exit = 0.0
        self.instThrust = 0.0

        # bookkeeping
        self.total_propellant = 0.0

    def inst(self, m_dot_ox, _):
        # -----------------------------
        # 1. Regression & fuel mass flow
        # -----------------------------
        G_ox = m_dot_ox / self.A_port
        self.r_dot_t = self.a * (G_ox**self.n)

        dV_fuel = np.pi * (
            (self.radius + self.r_dot_t * self.timestep) ** 2 - self.radius**2
        ) * self.L
        self.m_dot_fuel = self.rho_fuel * dV_fuel / self.timestep

        if self.m_fuel_t <= 0:
            self.m_dot_fuel = 0.0
            self.OF = 0.0
        else:
            self.OF = m_dot_ox / max(self.m_dot_fuel, 1e-9)

        self.m_dot_cc_t = m_dot_ox + self.m_dot_fuel

        # -----------------------------
        # 2. RocketCEA thermo calls
        # -----------------------------
        MW, gamma = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = R_UNIV / MW
        self.gamma = gamma

        h = self.C.get_Chamber_H(self.P_cc, self.OF, self.expratio)  # J/kg
        cp = self.C.get_Chamber_Cp(self.P_cc, self.OF, self.expratio, 0)  # J/kg-K
        self.cp = cp

        T_comb = self.C.get_Tcomb(self.P_cc, self.OF)

        # -----------------------------
        # 3. Energy balance (ΔU = in - out)
        # -----------------------------
        m_dot_in = self.m_dot_cc_t
        h_in = h
        m_dot_out = self.m_dot_cc_t
        h_out = h
        dU = (m_dot_in * h_in - m_dot_out * h_out) * self.timestep

        V_cc = self.L * self.A_port
        rho_cc = self.P_cc / (self.R * max(self.T_cc, 1.0))
        m_cc = rho_cc * V_cc

        dT = dU / (m_cc * cp) if m_cc > 0 else 0.0
        self.T_cc = self.T_cc + dT

        # -----------------------------
        # 4. Pc closure via choked-flow eqn
        # -----------------------------
        A_t = self.nozzle.A_throat
        self.P_cc = (
            self.m_dot_cc_t
            * np.sqrt(self.T_cc)
            / (
                A_t
                * np.sqrt(self.gamma / self.R)
                * (2.0 / (self.gamma + 1.0)) ** ((self.gamma + 1.0) / (2.0 * (self.gamma - 1.0)))
            )
        )

#BUG: THIS IS BAD
        if self.P_cc > 3.4e6:
            self.P_cc *=0.7

        # -----------------------------
        # 5. Nozzle expansion (always choked)
        # -----------------------------
        exit_mach = self.C.get_MachNumber(self.P_cc, self.OF, self.expratio, 0, 1)
        P_exit = self.P_cc / (
            (1 + (self.gamma - 1) / 2 * exit_mach**2) ** (self.gamma / (self.gamma - 1))
        )
        self.v_exit = np.sqrt(
            (2 * self.gamma)
            / (self.gamma - 1)
            * self.R
            * self.T_cc
            * (1 - (P_exit / self.P_cc) ** ((self.gamma - 1) / self.gamma))
        )

        # -----------------------------
        # 6. Thrust
        # -----------------------------
        self.instThrust = (self.m_dot_cc_t * self.v_exit) + self.nozzle.A_exit * (
            P_exit - self.P_atm
        )

        # -----------------------------
        # 7. Fuel regression update
        # -----------------------------
        self.m_fuel_t -= self.m_dot_fuel * self.timestep
        self.radius += self.r_dot_t * self.timestep
        self.A_port = np.pi * self.radius**2
        self.total_propellant += self.m_dot_cc_t * self.timestep

        return {"P_cc": self.P_cc, "thrust": self.instThrust}






"""
class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, L_star, m_fuel_i, rho_fuel, a, n, L, A_port, nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        self.C = C
        #initalizing y and R
        self.y = 0 
        self.R = 0 # J/Kg*K

        self.rho_fuel = rho_fuel 
        self.m_fuel_t = m_fuel_i 
        self.a = a 
        self.n = n 

        #fuel grain geometry
        self.L = L
        self.A_port = A_port

        self.expratio = self.nozzle.A_exit / self.nozzle.A_throat
        #print(self.expratio)

        self.timestep = timestep

        # initial values (treat as private members)
        self.P_cc = P_atm
        self.P_atm = P_atm
        self.v_exit = 0
        self.r_dot_t = 0.1 

        self.m_floss = 0
        self.m_dot_cc_t = 0        
        self.instThrust = 0

        self.OF = 0.5
        self.radius = np.sqrt(A_port / np.pi)

        # setup
        self.m_dot_fuel = 0
        self.total_propellant = 0

    def inst(self, m_dot_ox, _):
        # print( "ox m dot: ", m_dot_ox, " cc mass flow rate: ", self.m_dot_cc_t, " port area: ", self.A_port_t )
        # oxidizer mass flux iteration loop to get mass flow rate

        i = 0
        while i < 2:
            if i >= 1:
                G_ox_t = (m_dot_ox + self.m_dot_cc_t) / (2 * self.A_port)
            else:
                G_ox_t = m_dot_ox / self.A_port

            self.r_dot_t = self.a * (G_ox_t)**self.n
            
            #self.m_dot_fuel = self.rho_fuel * np.pi * ((self.r_dot_t + np.sqrt(self.A_port_t / np.pi))**2 - (self.A_port_t / np.pi)) * self.L
            #BUG: THE ABOVE IS WRONG, r_dot is bad!
            self.m_dot_fuel = self.rho_fuel * np.pi * ((self.r_dot_t*self.timestep + np.sqrt(self.A_port / np.pi))**2 - (self.A_port / np.pi)) * self.L
            
            if self.m_fuel_t <= 0:
                self.m_dot_fuel = 0
                self.OF = 0
            else:
                self.OF = m_dot_ox / self.m_dot_fuel

            i += 1

        self.m_dot_cc_t = m_dot_ox + self.m_dot_fuel

        # solve fluid properties
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K) #NOTE: IS THIS TOO BIG? UNIT ERROR RN?
        self.y = fluid_prop[1] # (-)

        #print(self.R, self.y)

        # call cea to get cc and exit temperature
        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.expratio, 0, 1)
        T_cc = temperatures[0]
        
        self.P_cc = (self.m_dot_cc_t /  self.nozzle.A_throat ) * np.sqrt( self.R*T_cc )  / ( np.sqrt( self.y * (2 / (self.y+1))**( (self.y+1)/(self.y-1) ) ) )
        print(self.P_cc,self.m_dot_cc_t,self.nozzle.A_throat,self.R,T_cc,self.y)


        #resolve fluid properties
        fluid_prop = self.C.get_Chamber_MolWt_gamma(self.P_cc, self.OF, self.expratio)
        self.R = 8314 / fluid_prop[0] # J/(kg K)
        self.y = fluid_prop[1] # (-)

        # call cea to get cc and exit temperature
        temperatures = self.C.get_Temperatures(self.P_cc, self.OF, self.expratio, 0, 1)
        T_cc = temperatures[0]

        #solve critical pressure:
        P_crit = self.P_cc*(1+ (self.y-1)/2)**((-1)*self.y/(self.y-1))


        
        #if subsonic at throat!
        if(P_crit < self.P_atm):
            #print("flow no longer choked! --> subsonic ")
            P_exit = self.P_atm
            T_exit = T_cc / (self.P_cc/P_exit)**((self.y-1)/self.y)
            #print(T_exit)
            cp = self.C.get_Chamber_Cp(self.P_cc, self.OF, self.expratio,0) #this is equilibrium check!
            self.v_exit = np.sqrt((2*cp)*(T_cc-T_exit))
            #print(self.v_exit, self.m_dot_cc_t)
        else:
        #supersonic #TODO: reverse case for speed? #TODO: make subsonic case a funciton where you pass in noteworthy upstream pressure
            if( 0.5 * self.P_atm < P_crit < 1.5* self.P_atm): #check that this is about equal and 50% deviation is about right
                #print("sonic throat subsonic exit case")
                P_exit = self.P_atm
                T_exit = T_cc / (P_crit/P_exit)**((self.y-1)/self.y) #pressure at throat will decrease from subsonic nozzle
                cp = self.C.get_Chamber_Cp(self.P_cc, self.OF, self.expratio,0) #this is equilibrium check!
                self.v_exit = np.sqrt((2*cp)*(T_cc-T_exit))

            else:
                exit_mach = self.C.get_MachNumber(self.P_cc, self.OF, self.expratio, 0, 1)

                P_exit = self.P_cc / (1 + ((self.y - 1) / 2) * (exit_mach**2) )**(self.y / (self.y - 1))

                self.v_exit = np.sqrt(((2 * self.y) / (self.y - 1)) * self.R * T_cc * (1 - (P_exit / self.P_cc)**((self.y - 1) / self.y)))
        



        
        self.instThrust = (self.m_dot_cc_t * self.v_exit) + self.nozzle.A_exit * (P_exit - self.P_atm)


        # Recalculate new fuel grain size
        self.m_fuel_t = self.m_fuel_t - self.m_dot_fuel * self.timestep

        self.total_propellant = self.total_propellant + self.m_dot_cc_t * self.timestep

        self.radius = self.radius + self.r_dot_t * self.timestep
        self.A_port = np.pi * self.radius**2

        return {"P_cc": self.P_cc, "thrust": self.instThrust}

"""

"""

class hybrid_cc_w_fuel_grain_model(BaseChamber):
    def __init__(self, L_star, m_fuel_i, rho_fuel, a, n, L, A_port, nozzle, P_atm, timestep, C):
        super().__init__(nozzle=nozzle, P_atm=P_atm, timestep=timestep, C=C)

        self.C = C

        # fuel grain properties
        self.rho_fuel = rho_fuel
        self.a = a
        self.n = n
        #self.h_latent_fuel = 400e3  # J/kg, placeholder

        # fuel grain geometry and mass
        self.m_fuel_t = m_fuel_i
        self.L = L
        #self.A_port = A_port
        self.r = np.sqrt(A_port / np.pi)

        print("r!", self.r)


        # oxidizer CEA object
        ox_card = """#fuel paraffin  N 2   O 1    wt%=100.00"""
"""        add_new_fuel('N2O', ox_card)
        self.oxCEA = CEA_Obj(oxName='N2O', fuelName='N2O')

        self.timestep = timestep
        self.m_dot_fuel = 0.0





    def inst(self, m_dot_ox, _):
        # --- 1. Regression
        A_port = np.pi * self.r**2
        G_ox = m_dot_ox / max(A_port, 1e-12)
        print("G_ox: ", G_ox)
        r_dot = self.a * G_ox**self.n
        m_dot_fuel = self.rho_fuel * (2*np.pi*self.r*self.L) * r_dot
        OF = m_dot_ox / max(m_dot_fuel, 1e-9)

        # --- 2. Residual closure
        def mass_balance(P_cc):
            T_cc = self.C.get_Tcomb(P_cc, OF)
            MW, gamma = self.C.get_Chamber_MolWt_gamma(P_cc, OF,
                                                       self.nozzle.expratio)
            R_spec = R_UNIV / MW

            #print("thermochemistry: ", P_cc, T_c, gamma, R_spec)
            _, m_dot_noz = self.nozzle.sol_thrust(P_cc, T_cc, gamma, R_spec)

            #print("res: ", ((m_dot_ox + m_dot_fuel) - m_dot_noz), m_dot_ox, m_dot_fuel, P_cc, OF)
            return (m_dot_ox + m_dot_fuel) - m_dot_noz

        sol = root_scalar(mass_balance,
                          bracket=[1e5, 1e8],
                          method="bisect")
        if not sol.converged:
            raise RuntimeError("Pcc solve failed")
        self.P_cc = sol.root

        # --- 3. Final thermo + thrust
        T_cc = self.C.get_Tcomb(self.P_cc, OF)
        MW, gamma = self.C.get_Chamber_MolWt_gamma(self.P_cc, OF,
                                                   self.nozzle.expratio)
        R_spec = R_UNIV / MW
        self.instThrust, _ = self.nozzle.sol_thrust(self.P_cc, T_cc, gamma, R_spec)

        # --- 4. Update port radius
        self.r += r_dot * self.timestep

        return {"P_cc": self.P_cc, "thrust": self.instThrust}
    
"""