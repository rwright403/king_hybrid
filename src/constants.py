import numpy as np

def intometer(x):
    return x*0.0254

###LAUNCH PAD DATA FOR LAUNCH CANADA
latitude=47.989083
longitude=-81.853361
elevation=370.3

launch_rail_length=25 #m
inclination = 85 #deg from ground
heading = 0

year = 2023
month = 10
date = 24
hour = 13

###ENGINE DATA
oxName = 'N2O'
rho_ox_liq = 1220 #kg/m^3
rho_ox_gas = 1.9277 #kg/m^3

fuelName = 'paraffin'
rho_fuel = 900 # kg/m^3

# RocketCEA doesnt have paraffin built in: CREATE IT BELOW
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

r_tank = intometer(5.5/2) #m --> for trajectory sim
height_tank = intometer(40) #m --> for trajectory sim
V_tank = 0.00937 #m^3

P_tank = 5000000 #Pa
fill_level = 0.6
C_inj = 0.00001
P_atm = 101325 #Pa

### Sim Variables
timestep = 0.05 #s
all_error = 0.01 
sim_time = 10 #s (time engine will be simulated over)

###ROCKET DATA --> MVH-1
rocket_fuselage_rad = intometer(5.5/2) #m --> for trajectory sim
rocket_dry_mass = 30 #kg

nosecone_shape = 'Power Series'
nosecone_length = 0.47 #m

###PRODUCE THRUST CURVE GRAPH
thrust_curve_graphs = True

###FILEPATHS FOR VALIDATION
model_thrust_file_path = r'./src/thrust.csv'
model_p_cc_file_path = r'./src/p_cc.csv'
model_p_tank_file_path = r'./src/p_tank.csv'

exp_thrust_file_path = r'./src/bens_validation_data/UofT_Deliverance_II/UofT_Deliverance_II_Thrust.csv'
exp_p_cc_file_path = r'./src/bens_validation_data/UofT_Deliverance_II/UofT_Deliverance_II_CC_Pressure.csv'
exp_p_tank_file_path = r'./src/bens_validation_data/UofT_Deliverance_II/UofT_Deliverance_II_Tank_Pressure.csv'

###SENSITIVITY ANALYSIS INFORMATION!!!! (only works for engine rn)
"""
Variables to analyize:
fill_level
C_inj
V_tank
P_tank
m_fuel_i
a
L
A_port_i
A_throat
A_exit
"""

test_var = "L"
min_bound = 0.2
max_bound = 0.4
num_iterations=5

#TODO: ADD n!!!!