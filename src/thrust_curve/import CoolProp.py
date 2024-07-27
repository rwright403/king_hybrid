import CoolProp.CoolProp as CP

# Upstream enthalpy of liquid
P_in = 5e6  # Pa
T_in = 284  # K
h_liquid = CP.PropsSI('H', 'P', P_in, 'T', T_in, 'N2O')  # Example: 300 kJ/kg = 300000 J/kg

#h_liquid is issue not using correct input pressure and temperature

# Downstream enthalpy of vapor
P_out = 1.5e6  # Pa
h_vapor_exit = CP.PropsSI('H', 'P', P_out, 'Q', 1, 'N2O')  # Example: 400 kJ/kg = 400000 J/kg

# Density of vapor at downstream pressure
rho_vapor_exit = CP.PropsSI('D', 'P', P_out, 'Q', 1, 'N2O')  # Example: 10 kg/m³

# Calculate mass flow rate
mass_flow_rate = rho_vapor_exit * 0.25 * 3.14 * 0.01 * 0.6 * (2 * (h_liquid - h_vapor_exit))**0.5

print(f'Mass Flow Rate: {mass_flow_rate} kg/s', rho_vapor_exit)