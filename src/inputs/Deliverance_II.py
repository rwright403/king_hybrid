"""Sim Variables"""

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


### CC models ###

""" 1 --> hybrid_cc_w_fuel_grain"""
#
oxName = 'N2O'
fuelName = 'paraffin'
timestep  = timestep
m_fuel_i = 1.5
rho_fuel = 900

# RocketCEA doesnt have paraffin built in: CREATE IT BELOW
#C32H66 from RPA Paraffin Wax Composition
CEA_fuel_str = f"""
#fuel paraffin  C 20   H 42    wt%=100.00
#h,KJ/mol=-13313.49  t(k)=298.15   rho,kg/m3={rho_fuel}
"""

a = 0.155/1000 #m/s
n = 0.5
L = 0.3852333 #m
A_port_i = 0.0038319753 #m^2
P_atm = P_atm
A_throat = 0.0010653525 #m^2
A_exit = 0.00531921243 #m^2

     
""" 2 --> adiabatic_lre_cc"""
#
"""
oxidizer_name = None,
fuel_name = None
A_throat = None
A_exit = None 
P_atm = P_atm
TIMESTEP = timestep
"""


### Tank models ###

""" 1 --> bens_ox_tank"""
#

oxName = 'N2O'
timestep = timestep
m_ox = 8.2 #kg
C_inj_1 = 0.0000294229787 #m^2
V_tank = 0.012 #m^3
P_tank = 3920207.656 #Pa
P_atm = P_atm 
all_error = all_error

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
"""
pressurant_name = None 
m_pressurant  = None 
fuel_name = None
m_fuel = None 
P_fueltank = None 
ID_PROPTANK = None 
TIMESTEP = timestep
"""


### Plotting
exp_thrust_file_path = None
exp_p_cc_file_paths = [None,None]
exp_p_tank_file_path = None