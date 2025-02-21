import matplotlib.pyplot as plt
import numpy as np
import csv
import importlib
import inspect

from src.models.bens_ox_tank import model

def read_csv(file_path):   
    x = []
    y = []
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:

            if not row:
                continue

            try:
                x.append(float(row[0]))
                y.append(float(row[1]))
            except ValueError:
                continue

    return np.array(x), np.array(y)

class cc_pressure():
    def __init__(self, exp_p_cc_file_path): #want tank temperature for coolprop
        self.Time_arr, self.Tank_Pres_arr = read_csv(exp_p_cc_file_path)

    def inst(self,t):
        #interpolate data from vectors to get real time pressure for spi model
        P_interpolated = np.interp(t, self.Time_arr, self.Tank_Pres_arr)
        return P_interpolated

def run_nitrous_validation(inputs):

    exp_time_p_ox_tank, exp_p_ox_tank = read_csv(inputs.exp_p_ox_tank_file_path)

    r1cc = cc_pressure(inputs.exp_p_cc_file_path)
    r1ox = model(inputs.oxName, inputs.TIMESTEP, inputs.m_ox, inputs.Cd_1, inputs.A_inj_1,
                inputs.V_tank, inputs.P_tank, inputs.P_atm, inputs.all_error, inputs.inj_model)
    
    t = 0
    TIMESTEP = inputs.TIMESTEP
    sim_time = inputs.sim_time

    time_arr = []
    p_ox_tank_arr = []
    m_arr = []
    m_dot_arr = []


    while t < sim_time:
        p = r1cc.inst(t)
        if p==0:
            p = 1e5
        r1ox.inst(p)

        time_arr.append(r1ox.t)
        p_ox_tank_arr.append(r1ox.P_tank)
        m_arr.append(r1ox.m_ox)
        m_dot_arr.append(r1ox.m_dot_ox)

        #print(r1ox.P_tank, p, t)
        #print(t, "\n")

        t += TIMESTEP


    plt.subplot(1,3,1)
    plt.axvline(x=4.1, color='red', linestyle='--', linewidth=2, label='liquid phase depletes')
    plt.plot(time_arr,p_ox_tank_arr, label = 'model out NOS Tank')
    plt.plot(exp_time_p_ox_tank, exp_p_ox_tank, label = 'experimental NOS Tank')
    
    exp_time_p_cc, exp_p_cc = read_csv(inputs.exp_p_cc_file_path)
    plt.plot(exp_time_p_cc, exp_p_cc, label = 'experimental NOS Tank Outlet Pressure')

    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (Pa)')
    plt.title('System Pressures Over Time')
    plt.grid(True)
    plt.legend()

    #mass vs time
    plt.subplot(1,3,2)
    plt.axvline(x=4.1, color='red', linestyle='--', linewidth=2, label='liquid phase depletes')
    plt.plot(time_arr,m_arr, label = 'model out NOS Tank')
    
    mass_path = r'./src/inputs/msc_test_cases/tomasz_test_case_m.csv'
    exp_time_mass, exp_mass = read_csv(mass_path)
    plt.plot(exp_time_mass, exp_mass, label = 'mass')
    plt.xlabel('Time (s)')
    plt.ylabel('Tank Mass (kg)')
    plt.title('Tank Mass Over Time')
    plt.grid(True)
    plt.legend()

    #mass flow rate vs time
    plt.subplot(1,3,3)
    plt.axvline(x=4.1, color='red', linestyle='--', linewidth=2, label='liquid phase depletes')
    plt.plot(time_arr,m_dot_arr, label = 'model out NOS Tank')
    
    mass_flow_rate_path = r'./src/inputs/msc_test_cases/tomasz_test_case_m_dot.csv'
    exp_time_p_cc, exp_p_cc = read_csv(mass_flow_rate_path)
    plt.plot(exp_time_p_cc, exp_p_cc, label = 'experimental NOS Tank')
    plt.xlabel('Time (s)')
    plt.ylabel('Mass Flow Rate (kg/s)')
    plt.title('Tank Mass Flow Rate')
    plt.grid(True)
    plt.legend()

    plt.show()


