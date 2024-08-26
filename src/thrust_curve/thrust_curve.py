from src import constants

from src.thrust_curve.bens_ox_tank import OxTank
from src.thrust_curve.combustion_chamber import cc
   
import matplotlib.pyplot as plt
import numpy as np
import csv

def to_csv(x_arr, y_arr, filename):
    combined_arr = list(zip(x_arr,y_arr))
    with open(r'./src/' + f'{filename}' + '.csv' , 'w', newline='') as file:
        #print(file)
        writer = csv.writer(file)
        writer.writerows(combined_arr)

#INPUT VARIABLES
time_arr = []
m_dot_arr = []
thrust_arr = []
p_cc_arr = []
p_tank_arr = []
P_cc = constants.P_atm

r1ox = OxTank(constants.oxName, constants.timestep, constants.m_ox, constants.C_inj,
               constants.V_tank, constants.P_tank, constants.P_atm, constants.all_error, constants.inj_model)


r1cc = cc(constants.oxName, constants.fuelName, constants.CEA_fuel_str, constants.m_fuel_i, 
          constants.rho_fuel, constants.a, constants.n, constants.L, constants.A_port_i, 
          constants.P_atm, constants.A_throat, constants.A_exit, constants.timestep,)

###ENTER THRUST CURVE
r1ox.inst(P_cc)
while (r1cc.OF != 0):
    #print(r1cc.OF)
    
    r1cc.inst(r1ox.m_dot_ox)
    r1ox.inst(r1cc.P_cc)

    #add to arrays
    time_arr.append(r1ox.t)
    m_dot_arr.append(r1ox.m_dot_ox)
    thrust_arr.append(r1cc.instThrust)
    p_cc_arr.append(r1cc.P_cc) #NOTE: convert to atmospheric depending on data
    p_tank_arr.append(r1ox.P_tank)

    total_propellant = r1cc.total_propellant

    print(r1cc.OF)

    #print(r1ox.P_tank, r1cc.P_cc,r1ox.P_tank- r1cc.P_cc,  )

    #print(r1ox.m_ox, r1cc.m_fuel_t)

    #print(r1ox.t, r1cc.v_exit,r1cc.m_dot_cc_t,r1cc.R)

    #print("time: ", r1ox.t)

print("\n", "### mass balance ###")
print("total propellant calculated through sim: ", total_propellant+r1ox.m_ox+r1cc.m_fuel_t, "total starting/input propellant: ", constants.m_ox+constants.m_fuel_i, "difference (conservation of mass): ",total_propellant+r1ox.m_ox+r1cc.m_fuel_t -constants.m_ox-constants.m_fuel_i)

###WRITE CSV FOR FLIGHT SIM AND VALIDATION
to_csv(time_arr,m_dot_arr, "m_dot_ox")
to_csv(time_arr,thrust_arr, "thrust")
to_csv(time_arr,p_cc_arr, "p_cc")
to_csv(time_arr,p_tank_arr, "p_tank")



###PLOTS
if constants.thrust_curve_graphs == True:
    plt.subplot(1,4,1)
    plt.plot(time_arr,m_dot_arr)
    plt.xlabel('Time (s)')
    plt.ylabel('m_dot_ox (kg/s)')
    plt.title('Mass Flow Rate Over Time')
    plt.grid(True)

    plt.subplot(1,4,2)
    plt.plot(time_arr,thrust_arr)
    plt.xlabel('Time (s)')
    plt.ylabel('Thrust (N)')
    plt.title('Thrust Curve')
    plt.grid(True)

    plt.subplot(1,4,3)
    plt.plot(time_arr,p_cc_arr)
    plt.xlabel('Time (s)')
    plt.ylabel('Chamber Pressure (Pa)')
    plt.title('Chamber Pressure Over Time')
    plt.grid(True)

    plt.subplot(1,4,4)
    plt.plot(time_arr,p_tank_arr)
    plt.xlabel('Time (s)')
    plt.ylabel('Tank Pressure (Pa)')
    plt.title('Tank Pressure Over Time')
    plt.grid(True)


    plt.show()

