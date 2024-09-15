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

oxidizer_name = 'N2O'
fuel_name = 'Ethanol'
A_throat = 0.00102028641 #m^2
A_exit = 0.00527334324 #m^2
P_atm = P_atm
TIMESTEP = timestep



### Tank models ###

""" 1 --> bens_ox_tank"""
#

oxName = 'N2O'
timestep = timestep 
m_ox = 4.48 #kg 
#NOTE: GUESSING Cd
C_inj =  1* 0.00001735222#0.8* 0.0000136284 #(num_orifices * Cd * orifice_diam) Note: guessing Cd of 0.6, NOTE: when it doesnt work this is why :)
V_tank = 6.4e-3 # - from report: "5.8L of nos in a 6.4L tank"
P_tank = 5.171e6 #Pa
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

pressurant_name = 'N2' 
m_pressurant  = 0.12 #NOTE: guessing for now, they gave a volume, should i change inputs to this model?
fuel_name = 'Ethanol' #NOTE: This might not work, assuming 100% when they used 95% as well
m_fuel = 1.12 #kg 
P_fueltank = 8.21e6 #Pa
ID_PROPTANK = 0.0254*5 #m 
TIMESTEP = timestep
#NOTE: no V_tank????



### Plotting
exp_thrust_file_path = None
exp_p_cc_file_paths = [None,None]
exp_p_tank_file_path = None

### TODO: Add rocket definition