from src import constants

from src.thrust_curve.bens_ox_tank import OxTank
from src.thrust_curve.combustion_chamber import cc
   
import matplotlib.pyplot as plt
import numpy as np
import csv

#INPUT VARIABLES
time_arr = []
m_dot_arr = []
thrust_arr = []
P_cc = constants.P_atm

r1ox = OxTank(constants.oxName, constants.timestep, constants.fill_level, constants.C_inj,
               constants.V_tank, constants.P_tank, constants.P_atm, constants.all_error)


r1cc = cc(constants.oxName, constants.fuelName, constants.CEA_fuel_str, constants.m_fuel_i, 
          constants.rho_fuel, constants.a, constants.n, constants.L, constants.A_port_i, 
          constants.P_atm, constants.A_throat, constants.A_exit, constants.timestep)

###ENTER THRUST CURVE
while r1ox.t < 7:
    r1ox.inst(P_cc)
    r1cc.inst(r1ox.m_dot_ox, P_cc)

    time_arr.append(r1ox.t)
    m_dot_arr.append(r1ox.m_dot_ox)
    thrust_arr.append(r1cc.instThrust)

###WRITE CSV FOR FLIGHT SIM
m_dot_combined_arr = list(zip(time_arr,m_dot_arr))
thrust_combined_arr = list(zip(time_arr, thrust_arr))

with open(r'./src/m_dot_ox.csv' , 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(m_dot_combined_arr)

with open(r'./src/thrust.csv' , 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(thrust_combined_arr)


"""
###PLOTS
plt.subplot(1,2,1)
plt.plot(time_arr,m_dot_arr)
plt.xlabel('Time (s)')
plt.ylabel('m_dot_ox (kg/s)')
plt.title('Mass Flow Rate Over Time')
plt.grid(True)

plt.subplot(1,2,2)
plt.plot(time_arr,thrust_arr)
plt.xlabel('Time')
plt.ylabel('Thrust')
plt.title('Thrust Curve')
plt.grid(True)


plt.show()
"""
