analysis_mode = [1,1]

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

### SELECT INJECTOR MODEL!
# 1 --> SPI 
# 2 --> HEM
# 3 --> Dyer
inj_model = 1


### Sim Variables
timestep = 0.05 #s
all_error = 0.01 
sim_time = 8.5#s #UNUSED: TODO: Delete


###Other constants
#k_tank_wall = 0.237 #kW/(m K)



###Phoenix 1a
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