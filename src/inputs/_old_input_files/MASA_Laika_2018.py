"""Sim Variables"""
#https://www.soundingrocket.org/uploads/9/0/6/4/9064598/37_project_report.pdf
# analysis mode (select models for sim):
#                               oxidizer tank model 
#                              v
#analysis mode for liquid = [A,B,C]
#                            ^   ^
#                            |    fuel tank model (if liquid rocket engine, otherwise vector should be 2 numbers long w hybrid cc)
#                            cc model
analysis_mode = [2,1,3]

timestep = 0.05 #s
sim_time = 3.4 #s (time engine will be simulated over)
#NOTE: if this is too big and you are simulating over a range the script will break
#is it really a leading edge simulation software if the ux is poor?

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
oxidizer_name = oxidizer_name
fuel_name = fuel_name
A_throat = 0.00102028641 #m^2
A_exit = 0.00527334324 #m^2
P_atm = P_atm
TIMESTEP = timestep



### Tank models ###

""" 1 --> bens_ox_tank"""
#
oxName = oxidizer_name
timestep = timestep 
m_ox = 4.48 #kg 
#NOTE: GUESSING Cd
Cd_1 = 0.66
A_inj_1 = 0.00007471705 #m^2
V_tank = 6.4e-3 # - from report: "5.8L of nos in a 6.4L tank"
P_tank = 5.171e6 #Pa
P_atm = P_atm 
all_error = 0.01

### SELECT INJECTOR MODEL!
# 1 --> SPI 
# 2 --> HEM
# 3 --> Dyer
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

"""simpleAdiabaticPressurizedTank"""
#
pressurant_name = pressurant_name 
m_pressurant  = 0.12 #NOTE: estimated for now based on volume they gave in report, should i change inputs to this model?
fuel_name = fuel_name #NOTE: This might not work, assuming 100% when they used 95% as well
m_fuel = 1.12 #kg 
P_fueltank = 4.82633e6 #Pa
ID_PROPTANK = 0.0254*5 #m 
V_tank_2 = 2.16e-3 #m^3
Cd_2 = 0.62
A_inj_2 = 0.0000136284 #m^2
T_amb = T_amb
TIMESTEP = timestep


### Plotting
exp_thrust_file_path = r'./src/inputs/liquid_validation_data/MASA_Laika/MASA_Laika_Thrust.csv'
exp_p_cc_file_path = r'./src/inputs/liquid_validation_data/MASA_Laika/MASA_Laika_CC_Pressure.csv'
exp_p_ox_tank_file_path = r'./src/inputs/liquid_validation_data/MASA_Laika/MASA_Laika_Ox_Tank_Pressure.csv'
exp_p_fuel_tank_file_path = r'./src/inputs/liquid_validation_data/MASA_Laika/MASA_Laika_Fuel_Tank_Pressure.csv'

### Sensitivity Analysis:
test_var_name = "P_tank"
min_bound = 4.5e6
max_bound = 6.0e6
num_iterations = 3


"""
test_var_name = "A_inj_1"
min_bound = 0.00007471705*.8
max_bound = 0.00007471705*.2
num_iterations = 3
"""