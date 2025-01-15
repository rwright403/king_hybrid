from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical
import matplotlib.pyplot as plt
import traceback


def thermo_span_wagner(rho, T, param):
    # Constants for N2O
    R = 8.3144598 / 44.0128 * 1000 # Gas constant (kJ/kg*K)
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


# Global Constants:
R_U = 8.31446 #J/(mol K)
T_REF = 298.15 #K
P_REF = 101325 #Pa

n2o = Chemical('nitrous oxide', T=T_REF)

MW = (n2o.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o.omega - 0.26992*n2o.omega**2

#NOTE: T_min = 298 K, T_max = 1400 K
def solve_cp_g_nist(T):
    A = 27.67988
    B = 51.14898
    C = -30.64454
    D = 6.847911
    E = -0.157906
    F = 71.24934
    G = 238.6164
    H = 82.04824
    t = T/1000

    cp = A +B*t + C*t**2 + D*t**3 + E/t**2
    return cp/MW #convert to J/kg K


###NOTE: Ziliac Thesis cp eqn for T_min = 182K and T_max = 310K
def solve_cp_ig_ziliac(T):
    A = 21.62
    B = 72.81
    C = -57.78
    D = 18.3
    E = 0

    t = T/1000
    cp = (A +B*t + C*(t**2) + D*(t**3) + E/(t**2))/MW #NOTE: no 1000 term?
    return cp #J/kg K

T_1 = 298.15 #K
P_1 = 1e5 #2e6 #just for testing

cp_ig_ziliac = solve_cp_ig_ziliac(T_1)
n2o_ig_t = Chemical('N2O', T=T_1) 


cp_ig_thermo = n2o_ig_t.Cpg #HeatCapacityGas.T_dependent_property_integral(T_REF, T_1)/MW #+ 38.677054950656256/MW # convert J/mol/K to J/kg/K #NOTE: see which correlations are used here
preos_g = PR(Tc=n2o_ig_t.Tc, Pc=n2o_ig_t.Pc, omega=n2o_ig_t.omega, T=T_1, P=P_1)
cp_thermo = preos_g.Cp_dep_g/MW + cp_ig_thermo

#print("cp J/(mol K): ", cp_ig_ziliac*MW, cp_ig_thermo*MW, "(ziliac,thermo)")
#print("real nitrous cp from thermo: ", preos_g.Cp_dep_g/MW, " + ",  cp_ig_thermo, " = ", cp_thermo)
#print("trusting ziliac ideal gas for reference: ",  preos_g.Cp_dep_g/MW, " + ",  cp_ig_ziliac, " = ", (preos_g.Cp_dep_g/MW + cp_ig_ziliac))
#print("nist cpg: ", solve_cp_g_nist(T_1) )
#print(cp_ig_thermo)
#print("calculator: ", 207937.69043709588/878.7683344539829)
#print("percent error ziliac (accepted) and thermo: ", np.abs(cp_ig_thermo-cp_ig_ziliac)/cp_ig_ziliac*100)





#NOTE: T_min = 298 K, T_max = 1400 K
#https://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Mask=1
#NOTE: This returns enthalpy relative to T_REF = 298.15 and P_REF = 1 atm.
def solve_h_g_nist(T):
    A = 27.67988
    B = 51.14898
    C = -30.64454
    D = 6.847911
    E = -0.157906
    F = 71.24934
    #G = 238.6164
    H = 82.04824
    t = T/1000

    h = A*t + (B*t**2)/2 + (C*t**3)/3 + (D*t**4)/4 - E/t + F - H
    return (h/MW)*1000 #convert to J/kg from kJ/mol


T_2 = 298.15
P_2 = 1e5

#print(f"nist polynomial enthalpy: ", solve_h_g_nist(T_2))

#enthalpy check
n2o_ig_h = Chemical('N2O', T=T_2, P=P_2) 
preos_g = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_2, P=P_2)

h_ig_thermo = n2o_ig_h.H/MW
#h_thermo = h_ig_thermo + preos_g.H_dep_g/MW #NOTE: i think I am mixing ideal gas with non ideal gas properties!
#print("thermo enthalpy: ", h_ig_thermo, ' + ', preos_g.H_dep_g/MW, ' = ', h_thermo)



#want to test the CEOS object:
from thermo import ChemicalConstantsPackage, PRMIX, CEOSLiquid, CEOSGas, FlashPureVLS
fluid = "N2O"

constants, correlations = ChemicalConstantsPackage.from_IDs([fluid])
eos_kwargs = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas)

constants, correlations = ChemicalConstantsPackage.from_IDs(["N2O"])

gas = CEOSGas(PRMIX, T=T_2, P=P_2, zs=[1], HeatCapacityGases=correlations.HeatCapacityGases,
              eos_kwargs=eos_kwargs)
#print("enthalpy of CEOSGas: ", gas.H()/MW)



# Define the substance
substance = 'NitrousOxide'

# Create a range of temperatures for the saturated region
T_min, T_max = CP.PropsSI('Tmin', substance), CP.PropsSI('Tcrit', substance)
temperatures = np.linspace(T_min, T_max, 500)

# Initialize arrays for saturated properties
P_sat = np.zeros_like(temperatures)
V_sat_liq = np.zeros_like(temperatures)
V_sat_vap = np.zeros_like(temperatures)
S_sat_liq = np.zeros_like(temperatures)
S_sat_vap = np.zeros_like(temperatures)

# Calculate saturated properties
for i, T in enumerate(temperatures):
    P_sat[i] = CP.PropsSI('P', 'T', T, 'Q', 0, substance)
    V_sat_liq[i] = 1 / CP.PropsSI('D', 'T', T, 'Q', 0, substance)
    V_sat_vap[i] = 1 / CP.PropsSI('D', 'T', T, 'Q', 1, substance)
    S_sat_liq[i] = CP.PropsSI('S', 'T', T, 'Q', 0, substance) / 1000  # Convert J/(kg*K) to kJ/(kg*K)
    S_sat_vap[i] = CP.PropsSI('S', 'T', T, 'Q', 1, substance) / 1000  # Convert J/(kg*K) to kJ/(kg*K)

# Create subplots for P-V and T-S diagrams side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))

# Plot P-V diagram
ax1.plot(V_sat_liq, P_sat, 'b-', label='Saturated Liquid')
ax1.plot(V_sat_vap, P_sat, 'r-', label='Saturated Vapor')
ax1.set_title('Pressure-Volume Diagram for N2O')
ax1.set_xlabel('Volume (m^3/kg)')
ax1.set_ylabel('Pressure (Pa)')
#ax1.set_yscale('log')

#

T_h = 500 #K
P_h = 3.5e5 #Pa
size = 30
"""
preos = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_h, P=P_h)
print(preos.phase)
#solving cp:
n2o_ig_ref = Chemical('N2O', T=T_REF)
cp_ig_thermo = n2o_ig_ref.Cpg 
preos = PR(Tc=n2o_ig_ref.Tc, Pc=n2o_ig_ref.Pc, omega=n2o_ig_ref.omega, T=T_REF, P=P_REF)
cp_thermo_ref = preos.Cp_dep_g/MW + cp_ig_thermo


n2o_ig_h = Chemical('N2O', T=T_h) 
cp_ig_thermo = n2o_ig_h.Cpg
preos = PR(Tc=n2o_ig_h.Tc, Pc=n2o_ig_h.Pc, omega=n2o_ig_h.omega, T=T_h, P=P_h)
cp_thermo = preos.Cp_dep_g/MW + cp_ig_thermo

h_test_3 = preos.H_dep_l/MW + cp_thermo*T_h - cp_thermo_ref*T_REF + 0
print("enthalpy? ", preos.H_dep_l/MW, " + (", cp_thermo*T_h, " - ",  cp_thermo_ref*T_REF, ") = ", h_test_3)
"""




###note: need to find a pure gas, a pure liquid state and test, then use that method on metastable



preos = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_h, P=P_h)
print(preos.phase)

#solving cp:
n2o_ig_ref = Chemical('N2O', T=T_REF) 
preos = PR(Tc=n2o_ig_ref.Tc, Pc=n2o_ig_ref.Pc, omega=n2o_ig_ref.omega, T=T_REF, P=P_REF)
cp_thermo_ref = preos.Cp_dep_g/MW + n2o_ig_ref.Cpg

n2o_ig_h = Chemical('N2O', T=T_h) 
preos = PR(Tc=n2o_ig_h.Tc, Pc=n2o_ig_h.Pc, omega=n2o_ig_h.omega, T=T_h, P=P_h)
cp_thermo = preos.Cp_dep_g/MW + n2o_ig_h.Cpg

#NOTE: this way makes the most sense to me?
h_test_3 = preos.H_dep_g/MW + cp_thermo*T_h - cp_thermo_ref*T_REF
print("enthalpy? ", preos.H_dep_g/MW, " + (", cp_thermo*T_h, " - ",  cp_thermo_ref*T_REF, ") = ", h_test_3)

preos_1 = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_h, P=P_h)
D_h = (1/(preos_1.V_g/MW))
print("density input to bens: ", D_h)
#NOTE: seems like enthalpy conversion ben used might be part of the problem?
print("bens span wagner implementation: ", thermo_span_wagner(D_h, T_h, "h") + 7.3397e+05 - 467599.224397) #Convert from Span-Wagner enthalpy convention to NIST
h_g_nist = solve_h_g_nist(T_h)
print("nist polynomial for ideal gas: ", h_g_nist, h_g_nist+preos.H_dep_g/MW)
print("checking departure function: ", preos.H_dep_g/MW)




eos_specific_volume = np.linspace( (1/1000), (1/300), size)
solved_p = np.array([])
preos_l_v_arr = np.array([])
v_arr = np.array([])
coolpresarr = np.array([])

for i in eos_specific_volume:
    
    
    preos_1 = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_h, V= i*MW)
    pressure = preos_1.P
    
    preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_h, P=pressure)
    preos_l_v_arr = np.append(preos_l_v_arr, preos_l.V_l/MW)

    v_arr = np.append(v_arr, i)

    try:
        coolpressure = CP.PropsSI('P', 'T', T_h, 'D', (1/i), "N2O")
        #print("Coolprop: ", CP.PhaseSI('T', T_h, 'D', (1/i), "N2O"), CP.PropsSI('H','T', T_1, 'D', (1/i), "N2O"))
    except:
        coolpressure = 0
    coolpresarr = np.append(coolpresarr, coolpressure)

    try:

        #solving cp:
        n2o_ig_h = Chemical('N2O', T=T_REF)
        cp_ig_thermo = n2o_ig_h.Cpg 
        preos_g = PR(Tc=n2o_ig_h.Tc, Pc=n2o_ig_h.Pc, omega=n2o_ig_h.omega, T=T_REF, P=P_REF)
        cp_thermo_ref = preos_g.Cp_dep_g/MW + cp_ig_thermo


        n2o_ig_h = Chemical('N2O', T=T_h) 
        cp_ig_thermo = n2o_ig_h.Cpg
        preos_g = PR(Tc=n2o_ig_h.Tc, Pc=n2o_ig_h.Pc, omega=n2o_ig_h.omega, T=T_h, P=P_h)
        cp_thermo = preos_g.Cp_dep_g/MW + cp_ig_thermo

        #testing enthalpy:
        h_test_1 = preos_l.H_dep_l

        #h_test_2 = #n2o_ig_h.HeatCapacityGas.T_dependent_property_integral(0, T_h)/MW # convert J/mol/K to J/kg/K Note: enthalpy is a T and P dependent property
        #print(h_test_1)
        
        h_test_2 = cp_thermo*T_h - cp_thermo_ref*T_REF + 0
        

        #n2o_ig_ref = Chemical('N2O', T=T_REF) 
        h_test_3 = preos_l.H_dep_l/MW + cp_thermo*T_h - cp_thermo_ref*T_REF + 0

        h_test_4 = CP.PropsSI('H', 'T', T_h, 'D', (1/i), "N2O")
        #print("h_test_4: ", h_test_4)

        h_test_5 = thermo_span_wagner((1/i), T_h, "h") + 7.3397e+05 #Convert from Span-Wagner enthalpy convention to NIST
        #print("h_test_5: ", h_test_5)

        #print("h tests 1,2,3,4,5: ",  h_test_1, h_test_2, h_test_3, h_test_4, h_test_5, CP.PhaseSI('T', T_h, 'D', (1/i), "N2O"), f"(pressure: {pressure})")
        #print("h 1 test: ", h_test_1, f"(pressure: {pressure})")

        #testing cp:
        #preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_f, P=P_tank)
        #Cp_2 = preos_l.Cp_dep_l/MW + n2o.Cpg #J/K
    except Exception as e:
        print(f"bad! {e}\n(pressure: {pressure})")
        pressure = 0

    solved_p = np.append( solved_p, pressure)


"""latent_heat_evap = preos_l.Hvap(T_liq)*MW

preos_l = PR(Tc=n2o.Tc, Pc=n2o.Pc, omega=n2o.omega, T=T_f, P=P_tank)
Cp_2 = preos_l.Cp_dep_l/MW + n2o.Cpg #J/K"""



#ax1.plot(eos_specific_volume, solved_p, label='Isotherm')
ax1.plot(preos_l_v_arr, solved_p, label='Thermo Isotherm')
ax1.plot(v_arr, coolpresarr, label='Coolprop Equillibrium')



ax1.set_xscale('log')
ax1.legend(loc='best')
ax1.grid(True, which="both", ls="--")




# Plot T-S diagram
ax2.plot(S_sat_liq, temperatures, 'b-', label='Saturated Liquid')
ax2.plot(S_sat_vap, temperatures, 'r-', label='Saturated Vapor')
ax2.set_title('Temperature-Entropy Diagram for N2O')
ax2.set_xlabel('Entropy (kJ/(kg*K))')
ax2.set_ylabel('Temperature (K)')
ax2.legend(loc='best')
ax2.grid(True, which="both", ls="--")

plt.tight_layout()
plt.show()
