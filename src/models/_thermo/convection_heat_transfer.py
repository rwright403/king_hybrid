import numpy as np
import CoolProp.CoolProp as CP
from thermo import Chemical # I had issues with Coolprop's N2O viscosity model
from src.models._thermo.n2o_viscosity_polynomials import *
from src.models._thermo.n2o_thermo_span_wagner_class import lightweight_span_wagner_eos_cp, lightweight_span_wagner_eos_d_rho_dT_P

g = 9.81 #m/s^2

#NOTE: SIGN CONVENTION: Q_dot dir: (+) T_1 --> T_2 (f is fluid)
def solve_Q_dot_natural_convection_liq(rho_f, T_1, T_2, T_f, P_f, c, n, L, Area, fluid):
    if fluid == "N2O":
        print("fluid", fluid)
        n2o = Chemical('N2O', T=T_f, P=P_f) 
        k_f = n2o.kl # Conductivity W/(m K)
        dyn_visc_f = liq_dynamic_visc_polynomial(T_f)
        visc_f = dyn_visc_f/rho_f # Kinematic viscosity (m^2/s)

        cp_f = lightweight_span_wagner_eos_cp(rho_f, T_f) #(preos_l.Cp_dep_l/MW + cp_ig) #J/K 

        d_rho_dT_P = lightweight_span_wagner_eos_d_rho_dT_P(rho_f, T_f) 
        beta = d_rho_dT_P/rho_f     #(1/rho_f)*d_rho_dT_P

    else: #USE COOLPROP OK
        k_f = CP.PropsSI('L', 'T', T_f, 'P', P_f, fluid)  # Conductivity W/(m K)
        dyn_visc_f = CP.PropsSI('V', 'T', T_f, 'P', P_f, fluid)  # Dynamic viscosity (Pa s)
        visc_f = dyn_visc_f/rho_f 
        
        cp_f = CP.PropsSI("Cpmass", "T", T_f, "P", P_f, fluid)
        beta = CP.PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T_f, "P", P_f, fluid)

    Gr = ((L**3)*g*beta*np.abs(T_2 - T_1) ) / (visc_f**2)
    Pr = (cp_f*visc_f)/ k_f
    X = Gr*Pr

    h = c * (k_f/L) * X**n

    Q_dot = h*Area*(T_1-T_2)

    return Q_dot #NOTE: Q_dot + going into (2)



#NOTE: SIGN CONVENTION: Q_dot dir: (+) T_1 --> T_2 
def solve_Q_dot_natural_convection_gas(rho_f, T_1, T_2, T_f, P_f, c, n, L, Area, fluid): #BUG: potential mistake, solving _1 properties with _2 inputs, double check this is likely a mistake
    if fluid == "N2O":
        n2o = Chemical('N2O', T=T_f, P=P_f)  #TODO: units here!!!
        k_f = n2o.kg # Conductivity W/(m K)
        dyn_visc_f = gas_dynamic_visc_polynomial(T_f)
        visc_f = dyn_visc_f/rho_f
        #visc_f = get_n2o_viscosity(T_f, P_f, "vapor") # Kinematic viscosity (m^2/s) # _old

        cp_f = lightweight_span_wagner_eos_cp(rho_f, T_f) #(preos_l.Cp_dep_l/MW + cp_ig) #J/K 

        d_rho_dT_P = lightweight_span_wagner_eos_d_rho_dT_P(rho_f, T_f) 
        beta = d_rho_dT_P/rho_f     #(1/rho_f)*d_rho_dT_P

    else: #USE COOLPROP OK
        k_f = CP.PropsSI('L', 'T', T_f, 'P', P_f, fluid)  # Conductivity W/(m K)
        dyn_visc_f = CP.PropsSI('V', 'T', T_f, 'P', P_f, fluid)  # Dynamic viscosity (Pa s)
        visc_f = dyn_visc_f/rho_f 
        
        cp_f = CP.PropsSI("Cpmass", "T", T_f, "P", P_f, fluid)
        beta = CP.PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T_f, "P", P_f, fluid)

    Gr = ((L**3)*g*beta*np.abs(T_2 - T_1) ) / (visc_f**2)
    Pr = (cp_f*visc_f)/ k_f
    X = Gr*Pr

    h = c * (k_f/L) * X**n

    Q_dot = h*(Area)*(T_1-T_2)

    return Q_dot