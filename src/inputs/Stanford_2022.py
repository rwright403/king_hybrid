import numpy as np

"""Sim Variables"""
#https://www.soundingrocket.org/uploads/9/0/6/4/9064598/37_project_report.pdf TODO: UPDATE
# analysis mode (select models for sim):
#                               oxidizer tank model 
#                              v
#analysis mode for liquid = [A,B,C]
#                            ^   ^
#                            |    fuel tank model (if liquid rocket engine, otherwise vector should be 2 numbers long w hybrid cc)
#                            cc model
analysis_mode = [2,1,3]

timestep = 0.05 #s
sim_time = 15 #s (time engine will be simulated over)
#NOTE: if this is too big and you are simulating over a range the script will break
#is it really a leading edge simulation software if the ux isnt poor?

### PROGRAM OUTPUT:
thrust_curve_graphs = True

### ENVIRONMENTAL DATA
P_atm = 101325 #Pa
T_amb = 293.15 #K

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
#fuel_name = Kerosene
pressurant_name = 'N2' 

### CC models ###

""" 1 --> hybrid_cc_w_fuel_grain"""
#"""fuck hybrids"""
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
oxidizer_name = oxidizer_name
fuel_name = 'Kerosene'
A_throat = 0.25*np.pi*((0.0254*1)**2) #m^2
A_exit = 3*A_throat #m^2
P_atm = P_atm
TIMESTEP = timestep



### Tank models ###

""" 1 --> bens_ox_tank"""
#
oxName = oxidizer_name
timestep = timestep 
m_ox = 10#kg 
Cd_1 = 1
A_inj_1 = 27.08e-6 #m^2
V_tank = (0.0254*42.47) *0.25*np.pi*((0.0254*4.75)**2) #m^3
P_tank = 6894.76* 590 #Pa
P_atm = P_atm 
all_error = 0.01

### SELECT INJECTOR MODEL!
# 1 --> SPI 
# 2 --> HEM
# 3 --> Dyer
# 4 --> just use 4
inj_model = 1


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

""" --> simpleAdiabaticPressurizedTank"""
#
pressurant_name = pressurant_name
m_pressurant  = 0.12 #NOTE: estimated for now based on amount of pressurant used by MASA's Laika
fuel_name_1 = "n-Dodecane"
#NOTE: NEED TO CHANGE THE S1 DELCARATION TO MAKE IT DIFFERENT SINCE COOLPROP AND ROCKETCEA DON'T HAVE SAME STRING AND SLIGHTLY DIFF FLUID FOR KEROSENE
m_fuel = 1.75 #kg 
P_fueltank = 6894.76* 610 #Pa
ID_PROPTANK = 0.0254*4.75 #m 
V_tank_2 = (0.0254*13.90) * 0.25*np.pi*((0.0254*4.75)**2) #m^3 
Cd_2 = 1
A_inj_2 = 2.4e-6 #m^2
T_amb = T_amb
TIMESTEP = timestep


### Calling experimental data for thrust curve
exp_thrust_file_path = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_Thrust.csv'
exp_p_cc_file_path = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_CC_Pressure.csv'
exp_p_ox_tank_file_path = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_Ox_Tank_Pressure.csv'
exp_p_fuel_tank_file_path = r'./src/inputs/liquid_validation_data/Stanford_2022/Stanford_2022_Fuel_Tank_Pressure.csv'


### Sensitivity Analysis:
test_var_name = "Cd_1"
min_bound = 0.5
max_bound = 0.7
num_iterations = 3


### TODO: Add rocket definition


### Estimates for Sizing Wizard ###
min_TW_ratio = 11
Cd_est = 0.6
mass_fraction_estimate = 0.2657741984
characteristic_len = 1 #m