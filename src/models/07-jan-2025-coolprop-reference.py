import CoolProp.CoolProp as CP

# Fluid
fluid = 'N2O'

# Triple point (CoolProp reference)
T_triple = 182.33  # K
P_triple = 0.00013 * 1e6  # Pa

# To shift the reference point, subtract h_new_ref from future calculations
def convert_enthalpy_to_nist_convention(T, P):
    T_ref = 298.15  # K
    P_ref = 101325 # Pa
    h_ref = CP.PropsSI('H', 'T', T_ref, 'P', P_ref, fluid)

    h = CP.PropsSI('H', 'T', T, 'P', P, fluid)
    return h - h_ref

# Example: Calculate enthalpy at T = 250 K and P = 1 atm with the adjusted reference
T_example = 600  # K
P_example = 101325 # Pa
h_adjusted = convert_enthalpy_to_nist_convention(T_example, P_example)
print(f"Adjusted enthalpy at T = {T_example} K and P = {P_example} Pa: {h_adjusted} J/kg")