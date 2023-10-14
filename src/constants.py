
oxName = 'N2O'
fuelName = 'paraffin'
rho_fuel = 900 # kg/m^3

# RocketCEA doesnt have paraffin built in: create it below
#C32H66 from RPA Paraffin Wax Composition
CEA_fuel_str = f"""
fuel paraffin  C 32   H 66    wt%=100.00
h,KJ/Kgmol=-1860600     t(k)=298.15   rho,kg/m3={rho_fuel}
"""

m_fuel_i = 0.8142857143 #kg
a = 0.155/1000 #m/s
n = 0.5 
L = 0.285 #m
A_port_i = 0.0009372051412 #m^2
 
A_throat = 0.0003433998524 #m^2
A_exit = 0.004229867915 #m^2

V_tank = 0.00937 #m^3
P_tank = 5000000 #Pa
fill_level = 0.6
C_inj = 0.0000083

P_atm = 101325 #Pa
timestep = 0.05 #s
all_error = 0.01 
