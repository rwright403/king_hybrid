from src import constants

from src.thrust_curve.bens_ox_tank import OxTank
from src.thrust_curve.combustion_chamber import cc
   
import matplotlib.pyplot as plt
import numpy as np


class DataClass:
    def __init__(self, time_arr, m_dot_arr, thrust_arr, p_cc_arr, p_tank_arr):
        self.time_arr_ = time_arr
        self.m_dot_arr_ = m_dot_arr
        self.thrust_arr_ = thrust_arr
        self.p_cc_arr_ = p_cc_arr
        self.p_tank_arr_ = p_tank_arr


#TODO: start by calculating variance

def sensitivityAnalysis(fill_level,C_inj,V_tank,P_tank,m_fuel_i, a, n, L, A_port_i, A_throat, A_exit):

    time_arr = []
    m_dot_arr = []
    thrust_arr = []
    p_cc_arr = []
    p_tank_arr = []

    P_cc = constants.P_atm

    r1ox = OxTank(constants.oxName, constants.timestep, fill_level, C_inj, V_tank, P_tank, constants.P_atm, constants.all_error)

    r1cc = cc(constants.oxName, constants.fuelName, constants.CEA_fuel_str, m_fuel_i, 
        constants.rho_fuel, a, n, L, A_port_i,constants.P_atm, A_throat, A_exit, constants.timestep)
        
    while r1ox.t < constants.sim_time:
        r1ox.inst(P_cc)
        r1cc.inst(r1ox.m_dot_ox, P_cc)

        time_arr.append(r1ox.t)
        m_dot_arr.append(r1ox.m_dot_ox)
        thrust_arr.append(r1cc.instThrust) #this is returning 0 every time?? its because self.inst thrust = 0!
        p_cc_arr.append(r1cc.P_cc)
        p_tank_arr.append(r1ox.P_tank)

        iteration = DataClass(time_arr,m_dot_arr,thrust_arr,p_cc_arr,p_tank_arr)

    return iteration

def produce_graphs(big_data,i_arr):
    j=1
    while(j<5):

        if(j==1):
            k=0
            for s in big_data:
                plt.subplot(1,4,j)
                plt.plot(s.time_arr_,s.m_dot_arr_, label= f'{i_arr[k]}')
                k+=1
            plt.legend()
            plt.xlabel('Time (s)')
            plt.ylabel('m_dot_ox (kg/s)')
            plt.title('Mass Flow Rate Over Time')
            plt.grid(True)

        if(j==2):
            k=0
            for s in big_data:
                plt.subplot(1,4,j)
                plt.plot(s.time_arr_,s.thrust_arr_, label= f'{i_arr[k]}')
                k+=1
            plt.legend()
            plt.xlabel('Time (s)')
            plt.ylabel('Thrust (N)')
            plt.title('Thrust Curve')
            plt.grid(True)

        if(j==3):
            k=0
            for s in big_data:
                plt.subplot(1,4,j)
                plt.plot(s.time_arr_,s.p_cc_arr_, label= f'{i_arr[k]}')
                k+=1
            plt.legend()
            plt.xlabel('Time (s)')
            plt.ylabel('Chamber Pressure (Pa)')
            plt.title('Chamber Pressure Over Time')
            plt.grid(True)

        if(j==4):
            k=0
            for s in big_data:
                plt.subplot(1,4,j)
                plt.plot(s.time_arr_,s.p_tank_arr_, label= f'{i_arr[k]}')
                k+=1
            plt.legend()
            plt.xlabel('Time (s)')
            plt.ylabel('Tank Pressure (Pa)')
            plt.title('Tank Pressure Over Time')
            plt.grid(True)
            

        j += 1


def update_i(i):
    i_arr.append(i)
    return i + (constants.max_bound-constants.min_bound)/(constants.num_iterations -1)


i = constants.min_bound
big_data = []
i_arr = []
#malding if statements, you hate to see it for each input, then call function with sensitivity analysis inside

if constants.test_var=="fill_level":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(i,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if constants.test_var=="C_inj":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.fill_level,i,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if constants.test_var=="V_tank":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.fill_level,constants.C_inj,i,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if constants.test_var=="P_tank":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,i,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if constants.test_var=="m_fuel_i":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,i,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if constants.test_var=="a":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,i,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if constants.test_var=="n":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,i,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if constants.test_var=="L":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,i,constants.A_port_i, constants.A_throat, constants.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if constants.test_var=="A_port_i":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,i, constants.A_throat, constants.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if constants.test_var=="A_throat":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, i, constants.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if constants.test_var=="A_exit":
    while(i<=constants.max_bound):
        big_data.append( sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, i) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)


#TODO: PROBABLY NICE TO PRINT OUT CONSTANTS
    
print("\nSensitivity Analysis Summary:\n")
print(f"test_var = \"{constants.test_var}\"")
print(f"min_bound = {constants.min_bound}")
print(f"max_bound = {constants.max_bound}")
print(f"num_iterations = {constants.num_iterations}\n")
print(f"ENGINE DATA\n")
print(f"oxName = '{constants.oxName}'")
print(f"rho_ox_liq = {constants.rho_ox_liq}")
print(f"rho_ox_gas = {constants.rho_ox_gas}")
print(f"fuelName = '{constants.fuelName}'")
print(f"rho_fuel = {constants.rho_fuel}")
print(f"m_fuel_i = {constants.m_fuel_i}")
print(f"a = {constants.a}")
print(f"n = {constants.n}")
print(f"L = {constants.L}")
print(f"A_port_i = {constants.A_port_i}")
print(f"A_throat = {constants.A_throat}")
print(f"A_exit = {constants.A_exit}")
print(f"r_tank = {constants.r_tank}")
print(f"height_tank = {constants.height_tank}")
print(f"V_tank = {constants.V_tank}")
print(f"P_tank = {constants.P_tank}")
print(f"fill_level = {constants.fill_level}")
print(f"C_inj = {constants.C_inj}")
print(f"P_atm = {constants.P_atm}")
print(f"timestep = {constants.timestep}")
print(f"sim_time = {constants.sim_time}")
#i love good code :)

plt.show()