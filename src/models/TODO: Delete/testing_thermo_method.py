from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical

# Global Constants:
T_REF = 298.15 #K
P_REF = 101325 #Pa

R_U = 8.31446 #J/(mol K) 

n2o = Chemical('nitrous oxide')

MW = (n2o.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o.omega - 0.26992*n2o.omega**2

TANK_DIAM = 0.0254*5.5 #m
CS_AREA = 0.25*np.pi*(TANK_DIAM**2) #m^2
g = 9.81 #m/s^2

TIMESTEP = 1e-4 #S



def solve_thermo_properties(T_liq, T_gas, m_liq, m_gas, m_dot_liq, m_dot_gas, V_TANK, V_liq_prev, V_gas_prev, rho_liq_prev, rho_gas_prev):

    # Initial guess for dV/dt_liq
    V_dot_liq_guess = 0.005  # You can set this based on the previous timestep or system behavior
    V_dot_liq_iter = 0.005
    
    volume_tolerance = 1e-6 
    pressure_tolerance = 1e4
    max_iterations = 200 #TODO: define on max V_dot? can we thermodynamically figure out a maximum with some sort of hand calc assumption?

    for i in range(max_iterations):
        # Calculate dV/dt_vap based on dV/dt_liq
        V_dot_gas = -V_dot_liq_guess

        #use initial guess to solve future V_liq and V_gas with fwd euler method:
        V_liq_est = V_liq_prev + TIMESTEP * V_dot_liq_guess
        V_gas_est = V_gas_prev + TIMESTEP * V_dot_gas
        
        # Calculate d_rho/dt for liquid and vapor
        rho_dot_liq = (1 / V_liq_est) * m_dot_liq - (m_liq / V_liq_est**2) * V_dot_liq_guess 
        rho_dot_gas = (1 / V_gas_est) * m_dot_gas - (m_gas / V_gas_est**2) * V_dot_gas       
        
        # Update densities using Forward Euler method
        rho_liq = rho_liq_prev + TIMESTEP * rho_dot_liq
        rho_gas = rho_gas_prev + TIMESTEP * rho_dot_gas
        
        # Calculate pressures using Peng-Robinson EOS #NOTE: try this first, if it doesnt look good try coolprop but that might not work not sure we can see 
        pr_eos_liq = PR(T=T_liq, V=(MW/rho_liq), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
        P_liq = pr_eos_liq.P
        pr_eos_gas = PR(T=T_gas, V=(MW/rho_gas), Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega)
        P_gas = pr_eos_gas.P
    
        
        # Check constraints:
        volume_constraint = abs(V_TANK - (m_liq / rho_liq + m_gas / rho_gas))
        pressure_constraint = abs(P_liq - P_gas)
        #TODO: pressure derivative constraint?
        
        # If both constraints are satisfied, return the solution
        if volume_constraint < volume_tolerance and pressure_constraint < pressure_tolerance:
            return P_liq, P_gas, rho_liq, rho_gas, V_dot_liq_guess

        V_dot_liq_guess += V_dot_liq_iter
        i+=1

    raise Exception("Failed to converge")