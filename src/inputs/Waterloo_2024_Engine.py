"""Sim Variables"""
import numpy as np

# analysis mode (select models for sim):
#                               oxidizer tank model 
#                              v
#analysis mode for liquid = [A,B,C]
#                            ^   ^
#                            |    fuel tank model (if liquid rocket engine, otherwise vector should be 2 numbers long w hybrid cc)
#                            cc model
analysis_mode = [1,1]

timestep = 0.05 #s
all_error = 0.01

### PROGRAM OUTPUT:
thrust_curve_graphs = True

### ENVIRONMENTAL DATA
P_atm = 101325 #Pa
T_amb = 275 #K

### Launch Canada Timmins Pad
latitude = 47.989083
longitude = -81.853361
elevation = 370.3

launch_rail_length =25 #m
inclination = 85 #deg from ground
heading = 0

year = 2023
month = 10
date = 24
hour = 13

### Propellant and Pressurants 
oxidizer_name = 'N2O'
fuel_name = 'Ethanol'
pressurant_name = 'N2' 

### CC models ###

""" 1 --> hybrid_cc_w_fuel_grain"""
#
"""
oxName = None
fuelName = None
timestep  = None
m_fuel_i = None
rho_fuel = None

# RocketCEA doesnt have paraffin built in: CREATE IT BELOW
#C32H66 from RPA Paraffin Wax Composition
#CEA_fuel_str = f"""
#fuel paraffin  C 20   H 42    wt%=100.00
#h,KJ/mol=-13313.49  t(k)=298.15   rho,kg/m3={rho_fuel}
"""

a = None
n = None
L = None
A_port_i = None
P_atm = P_atm
A_throat = None
A_exit = None
"""
     
""" 2 --> adiabatic_lre_cc"""
#

oxidizer_name = oxidizer_name,
fuel_name = fuel_name
A_throat = 0.25*np.pi*(0.0254*(1.975))**2
A_exit = A_throat * 4.1
P_atm = P_atm
TIMESTEP = timestep



### Tank models ###

""" 1 --> bens_ox_tank"""
#

oxName = oxidizer_name
timestep = timestep 
m_ox = None 
Cd_1 = 0.66 #guess, could be v wrong
A_inj_1 = 1 #m^2
V_tank = 0.25*np.pi*( (0.0254*5.625)**2 - (0.0254*2.84)**2 ) * (0.0254*64)
P_tank = 5.51581e6
P_atm = P_atm 
all_error = all_error 

### SELECT INJECTOR MODEL!
# 1 --> SPI 
# 2 --> HEM
# 3 --> Dyer
inj_model = 4


"""awful liquid"""
#
"""
propellant = None 
pressurant = None 
id_PROPTANK = None
P_proptank = None
m_prop = None 
V_PROPTANK = None 
TIMESTEP = timestep
pressurant_name = None 
P_prestank = None 
m_pres = None 
P_oxtank = None
V_PRESTANK = None
OUTLET_DIAM = None
"""

"""simpleAdiabaticPressurizedTank"""
#
"""
pressurant_name = pressurant_name 
m_pressurant  = 
fuel_name = fuel_name #NOTE: This might not work, assuming 100% when they used 95% as well
m_fuel = 4.25 #kg 
P_fueltank = 5.51581e6 #Pa
ID_PROPTANK = 0.0254*2.84 #m 
V_tank_2 = 0.25*np.pi*(0.0254*2.84**2 ) * (0.0254*64) #m^3
Cd_2 = 
A_inj_2 =  #m^2
T_amb = T_amb
TIMESTEP = timestep
"""

"""piping_real_fuel_tank_data"""
#

T_tank = T_amb
Cd_spi = 0.66 #guess, could be v wrong
print("NOTE*** COAXIAL INJECTOR UNKNOWN DISCHARGE COEFFS")
A_inj = None
fuel_name = fuel_name
fuel_tank_pressure_filepath = r'./src/inputs/liquid_validation_data/Waterloo_2024/Waterloo_2024_Fuel_Tank_Pressure.csv'




### Plotting
exp_thrust_file_path = r'./src/inputs/liquid_validation_data/MASA_Laika/MASA_Laika_Thrust.csv'
print("thrust not validated rn will be wrong")
exp_p_cc_file_path = r'./src/inputs/liquid_validation_data/Waterloo_2024/Waterloo_2024_CC_Pressure.csv'
exp_p_ox_tank_file_path = r'./src/inputs/liquid_validation_data/Waterloo_2024/Waterloo_2024_Oxidizer_Tank_Pressure.csv'
exp_p_fuel_tank_file_path = r'./src/inputs/liquid_validation_data/Waterloo_2024/Waterloo_2024_Fuel_Tank_Pressure.csv'

### TODO: Add rocket definition