import numpy as np
import CoolProp.CoolProp as CP
#import time
from thermo import Chemical # I had issues with Coolprop's N2O viscosity model
from src.models._thermo.n2o_viscosity_polynomials import *
from src.models._thermo.n2o_thermo_span_wagner_class import lightweight_span_wagner_eos_cp, lightweight_span_wagner_eos_d_rho_dT_P
from src.models._thermo.conductivity_lookup_tables.n2o_kg_lookup import N2OGasConductivityTable
kg = N2OGasConductivityTable()
from src.models._thermo.conductivity_lookup_tables.n2o_kl_lookup import N2OLiqConductivityTable
kl = N2OLiqConductivityTable()
from src.models._thermo.cp_lookup_tables.n2o_cp_gas_lookup import N2OCPGasTable
cp_gas_table = N2OCPGasTable()
from src.models._thermo.cp_lookup_tables.n2o_cp_liq_lookup import N2OCPLiqTable
cp_liq_table = N2OCPLiqTable()
from src.models._thermo.d_rho_dT_const_P_lookup_tables.n2o_d_rho_dT_const_P_gas_lookup import N2ODRHODTPGasTable
d_rho_dT_const_P_gas_table = N2ODRHODTPGasTable()
from src.models._thermo.d_rho_dT_const_P_lookup_tables.n2o_d_rho_dT_const_P_liq_lookup import N2ODRHODTPLiqTable
d_rho_dT_const_P_liq_table = N2ODRHODTPLiqTable()

g = 9.81 #m/s^2

#precompute air: USE COOLPROP OK
P_atm = 101325 #Pa
T_atm = 273.15+25 # K

k_f_air = CP.PropsSI('L', 'T', T_atm, 'P', P_atm, "Air")  # Conductivity W/(m K)
dyn_visc_f_air = CP.PropsSI('V', 'T', T_atm, 'P', P_atm, "Air")  # Dynamic viscosity (Pa s)
rho_f_air = CP.PropsSI('D', 'P', P_atm, 'T', T_atm, "Air")
visc_f_air = dyn_visc_f_air/rho_f_air 

cp_f_air = CP.PropsSI("Cpmass", "T", T_atm, "P", P_atm, "Air")
beta_air = CP.PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T_atm, "P", P_atm, "Air")


#NOTE: SIGN CONVENTION: Q_dot dir: (+) T_1 --> T_2 (f is fluid)
def solve_Q_dot_natural_convection_liq(rho_f, T_1, T_2, T_f, P_f, c, n, L, Area, fluid):
    if fluid == "N2O":
        #n2o = Chemical('N2O', T=T_f, P=P_f) 
        #k_f = n2o.kl # Conductivity W/(m K)
        k_f = kl.lookup(T=T_f, P=P_f)
        dyn_visc_f = liq_dynamic_visc_polynomial(T_f)
        visc_f = dyn_visc_f/rho_f # Kinematic viscosity (m^2/s)

        cp_f = cp_liq_table.lookup(rho_f, T_f) #lightweight_span_wagner_eos_cp(rho_f, T_f) #(preos_l.Cp_dep_l/MW + cp_ig) #J/K 

        d_rho_dT_P = d_rho_dT_const_P_liq_table.lookup(rho_f, T_f) #
        beta = d_rho_dT_P/rho_f     #(1/rho_f)*d_rho_dT_P

    else: #USE COOLPROP OK
        k_f = k_f_air #CP.PropsSI('L', 'T', T_f, 'P', P_f, fluid)  # Conductivity W/(m K)
        #dyn_visc_f = CP.PropsSI('V', 'T', T_f, 'P', P_f, fluid)  # Dynamic viscosity (Pa s)
        visc_f = visc_f_air #dyn_visc_f/rho_f 
        
        cp_f = cp_f_air #CP.PropsSI("Cpmass", "T", T_f, "P", P_f, fluid)
        beta = beta_air #CP.PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T_f, "P", P_f, fluid)

    Gr = ((L**3)*g*beta*np.abs(T_2 - T_1) ) / (visc_f**2)
    Pr = (cp_f*visc_f)/ k_f
    X = Gr*Pr

    h = c * (k_f/L) * X**n

    Q_dot = h*Area*(T_1-T_2)

    return Q_dot #NOTE: Q_dot + going into (2)



#NOTE: SIGN CONVENTION: Q_dot dir: (+) T_1 --> T_2 
def solve_Q_dot_natural_convection_gas(rho_f, T_1, T_2, T_f, P_f, c, n, L, Area, fluid): #BUG: potential mistake, solving _1 properties with _2 inputs, double check this is likely a mistake
    if fluid == "N2O":

        """
        n2o = Chemical('N2O', T=T_f, P=P_f)  #TODO: units here!!!
        k_f = n2o.kg # Conductivity W/(m K)
        """
        k_f = kg.lookup(T=T_f, P=P_f)

        dyn_visc_f = gas_dynamic_visc_polynomial(T_f)
        
        visc_f = dyn_visc_f/rho_f
        #visc_f = get_n2o_viscosity(T_f, P_f, "vapor") # Kinematic viscosity (m^2/s) # _old

        cp_f = cp_gas_table.lookup(rho_f, T_f) #lightweight_span_wagner_eos_cp(rho_f, T_f) #(preos_l.Cp_dep_l/MW + cp_ig) #J/K 

        d_rho_dT_P = d_rho_dT_const_P_gas_table.lookup(rho_f, T_f)
        beta = d_rho_dT_P/rho_f     #(1/rho_f)*d_rho_dT_P

    else: #USE COOLPROP OK
        k_f = k_f_air #CP.PropsSI('L', 'T', T_f, 'P', P_f, fluid)  # Conductivity W/(m K)
        #dyn_visc_f = CP.PropsSI('V', 'T', T_f, 'P', P_f, fluid)  # Dynamic viscosity (Pa s)
        visc_f = visc_f_air #dyn_visc_f/rho_f 
        
        cp_f = cp_f_air #CP.PropsSI("Cpmass", "T", T_f, "P", P_f, fluid)
        beta = beta_air #CP.PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T_f, "P", P_f, fluid)

    Gr = ((L**3)*g*beta*np.abs(T_2 - T_1) ) / (visc_f**2)
    Pr = (cp_f*visc_f)/ k_f
    X = Gr*Pr

    h = c * (k_f/L) * X**n

    Q_dot = h*(Area)*(T_1-T_2)

    return Q_dot