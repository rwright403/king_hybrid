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

fuelName = 'paraffin'
rho_fuel = 900 # kg/m^3

# RocketCEA doesnt have paraffin built in: CREATE IT BELOW
#C32H66 from RPA Paraffin Wax Composition
CEA_fuel_str = f"""
fuel paraffin  C 20   H 42    wt%=100.00
h,KJ/mol=-13313.49  t(k)=298.15   rho,kg/m3={rho_fuel}
"""

#heat of combustion -13313.49 found from:
# https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Resources/References/paraffin%20heat%20of%20formation.pdf
#

a = 0.155/1000 #m/s
n = 0.5

P_atm = 101325 #Pa
"""
m_fuel_i = 6.173983 #kg
L = 0.4 #m
A_port_i = 0.001963495408 #m^2
 
A_throat = 0.000697464985 #m^2
A_exit = 0.00417781528 #m^2

V_tank = 0.043 #m^3

P_tank = 6.5e+6 #Pa
m_ox = 18.30872404 #kg
"""

"""
num_orifices
Cd
Orifice_diam

"""
#C_inj = 0.00002243097155 #m^2 #(num_orifices * Cd * orifice_diam)


### SELECT INJECTOR MODEL!
# 1 --> SPI 
# 2 --> HEM
# 3 --> Dyer
inj_model = 1


### Sim Variables
timestep = 0.05 #s
all_error = 0.01 
sim_time = 8.5#s (time engine will be simulated over)
#NOTE: if this is too big and you are simulating over a range the script will break
#is it really a leading edge simulation software if it isnt annoying to use?

###Other constants
k_tank_wall = 0.237 #kW/(m K)


###ROCKET DATA FOR TRAJECTORY SIM--> MVH-1
rocket_fuselage_rad = intometer(5.5/2) #m --> for trajectory sim
rocket_dry_mass = 30 #kg

r_tank = intometer(5.5/2) #m --> for trajectory sim
height_tank = intometer(40) #m --> for trajectory sim

nosecone_shape = 'Power Series'
nosecone_length = 0.47 #m

###PRODUCE THRUST CURVE GRAPH
thrust_curve_graphs = True

###FILEPATHS FOR VALIDATION
model_thrust_file_path = r'./src/thrust.csv'
model_p_cc_file_path = r'./src/p_cc.csv'
model_p_tank_file_path = r'./src/p_tank.csv'

###Phoenix 1a
"""
name = "Phoenix 1a"

m_fuel_i = 6.173983 #kg
L = 0.4 #m
A_port_i = 0.001963495408 #m^2
 
A_throat = 0.000697464985 #m^2
A_exit = 0.00417781528 #m^2

V_tank = 0.043 #m^3

P_tank = 6.5e+6 #Pa
m_ox = 18.30872404 #kg
C_inj = 0.00002243097155 #m^2 (number of orifices * Cd * Orifice Area)


exp_thrust_file_path = r'./src/bens_validation_data/Phoenix_1a/Phoenix_1a_Thrust.csv'
exp_p_cc_file_path = r'./src/bens_validation_data/Phoenix_1a/Phoenix_1a_CC_Pressure.csv'
exp_p_tank_file_path = r'./src/bens_validation_data/Phoenix_1a/Phoenix_1a_Tank_Pressure.csv'
"""


###Deliverance II
name = "Deliverance II"

m_fuel_i = 1.5 #kg
L = 0.3852333 #m
A_port_i = 0.0038319753 #m^2

A_throat = 0.0010653525 #m^2
A_exit = 0.00531921243 #m^2

V_tank = 0.012 #m^3

P_tank = 3920207.656 #Pa
m_ox = 8.2 #kg
C_inj = 0.0000294229787 #m^2

exp_thrust_file_path = r'./src/bens_validation_data/UofT_Deliverance_II/UofT_Deliverance_II_Thrust.csv'
exp_p_cc_file_path = r'./src/bens_validation_data/UofT_Deliverance_II/UofT_Deliverance_II_CC_Pressure.csv'
exp_p_tank_file_path = r'./src/bens_validation_data/UofT_Deliverance_II/UofT_Deliverance_II_Tank_Pressure.csv'


###Boundless
"""
name = "Boundless"

m_fuel_i = 5.44311 #kg
L = 0.41130982 #m
A_port_i = 0.001302041 #m^2
 
A_throat = 0.001302041 #m^2
A_exit = 0.0061470149 #m^2

V_tank = 0.024956005 #m^3

P_tank = 4933902 #4.206e6 #Pa
m_ox = 14.61751413 #kg
C_inj = 0.00002364094414 #m^2

exp_thrust_file_path = r'./src/bens_validation_data/UofW_Boundless/UofW_Boundless_Thrust.csv'
exp_p_cc_file_path = r'./src/bens_validation_data/UofW_Boundless/UofW_Boundless_CC_Pressure.csv'
exp_p_tank_file_path = r'./src/bens_validation_data/UofW_Boundless/UofW_Boundless_Tank_Pressure.csv'
"""

### WATERLOO KISMET

###SENSITIVITY ANALYSIS INFORMATION!!!! (only works for engine rn)
"""
Variables to analyize:
fill_level TODO: update in sensitivity analysis
C_inj
V_tank
P_tank
m_fuel_i
a
n
L
A_port_i
A_throat
A_exit
"""

test_var = "C_inj"
min_bound = 1.71e-5
max_bound = 2.39e-5
num_iterations=3


