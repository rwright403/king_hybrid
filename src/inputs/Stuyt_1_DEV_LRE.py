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
fuel_name = fuel_name
A_throat = 0.001342491 #m^2
A_exit = 0.0055046156 #m^2
P_atm = P_atm
TIMESTEP = timestep



### Tank models ###

""" 1 --> bens_ox_tank"""
#

oxName = oxidizer_name
timestep = timestep 
m_ox = 7.650873122104811 #kg 
#NOTE: GUESSING Cd
C_inj_1 =  0.6 * 03.090605599220321e-5 #Note: guessing Cd of 0.6
V_tank = 0.01177057403 #m^3
P_tank = 5.2e6 #Pa
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
m_pressurant  = 0.12 #NOTE: estimated for now based on amount of pressurant used by MASA's Laika
fuel_name = fuel_name
m_fuel = 1.5301746244209624 #kg 
P_fueltank = 5.2e6 #Pa
ID_PROPTANK = 0.0254*5 #m 
V_tank_2 = 0.0037961342 #m^3 
C_inj_2 = 0.6*6.566075013100621e-6 #m^2
T_amb = T_amb
TIMESTEP = timestep


### Calling experimental data for thrust curve
exp_thrust_file_path = None
exp_p_cc_file_path = None
exp_p_ox_tank_file_path = None
exp_p_fuel_tank_file_path = None


### Sensitivity Analysis:
test_var_name = "C_inj_2"
min_bound = 0.5*6.566075013100621e-6 #m^2
max_bound = 0.7*6.566075013100621e-6 #m^2
num_iterations = 3


### TODO: Add rocket definition


### Estimates for Sizing Wizard ###
min_TW_ratio = 11
Cd_est = 0.6
mass_fraction_estimate = 0.2657741984
characteristic_len = 1 #m