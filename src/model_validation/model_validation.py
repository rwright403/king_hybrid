import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from src import inputs

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

def validate(inputs):

    #TODO:return summary of rocket!!!!
    #print("\n")
    #print("Rocket: ", inputs.name)

    if len(inputs.analysis_mode) == 2:
        
        ###READ FROM FILES!!!!
        model_time_thrust, model_thrust = read_csv(inputs.model_thrust_file_path)
        model_time_p_cc, model_p_cc = read_csv(inputs.model_p_cc_file_path)
        model_time_p_tank, model_p_tank = read_csv(inputs.model_p_tank_file_path)

        exp_time_thrust, exp_thrust = read_csv(inputs.exp_thrust_file_path)
        exp_time_p_cc, exp_p_cc = read_csv(inputs.exp_p_cc_file_path)
        exp_time_p_tank, exp_p_tank = read_csv(inputs.exp_p_tank_file_path)

        ###PLOT!!!
        plt.subplot(1,3,1)
        plt.plot(model_time_thrust, model_thrust, label='model output')
        plt.plot(exp_time_thrust, exp_thrust,label='data') 
        plt.xlabel('Time (s)')
        plt.ylabel('Thrust (N)')
        plt.title('Thrust Curve Validation')
        plt.grid(True)

        plt.subplot(1,3,2)
        plt.plot(model_time_p_cc, model_p_cc, label='model output')
        plt.plot(exp_time_p_cc, exp_p_cc, label='data')
        plt.xlabel('Time (s)')
        plt.ylabel('Chamber Pressure (Pa)')
        plt.title('Chamber Pressure Validation')
        plt.grid(True)

        plt.subplot(1,3,3)
        plt.plot(model_time_p_tank, model_p_tank, label='model output')
        plt.plot(exp_time_p_tank, exp_p_tank, label='data')
        plt.xlabel('Time (s)')
        plt.ylabel('Tank Pressure (Pa)')
        plt.title('Tank Pressure Validation')
        plt.grid(True)

        plt.legend()
        plt.show()

    if len(inputs.analysis_mode) == 3:
        
        ###READ FROM FILES!!!!
        model_time_thrust, model_thrust = read_csv(r'./src/thrust.csv')
        model_time_p_cc, model_p_cc = read_csv(r'./src/p_cc.csv')
        model_time_p_ox_tank, model_p_ox_tank = read_csv(r'./src/p_ox_tank.csv')
        model_time_p_fuel_tank, model_p_fuel_tank = read_csv(r'./src/p_fuel_tank.csv')

        exp_time_thrust, exp_thrust = read_csv(inputs.exp_thrust_file_path)
        exp_time_p_cc, exp_p_cc = read_csv(inputs.exp_p_cc_file_path)
        exp_time_p_ox_tank, exp_p_ox_tank = read_csv(inputs.exp_p_ox_tank_file_path)
        exp_time_p_fuel_tank, exp_p_fuel_tank = read_csv(inputs.exp_p_fuel_tank_file_path)

        ###PLOT!!!
        plt.subplot(1,4,1)
        plt.plot(model_time_thrust, model_thrust, label='model output')
        plt.plot(exp_time_thrust, exp_thrust,label='data') 
        plt.xlabel('Time (s)')
        plt.ylabel('Thrust (N)')
        plt.title('Thrust Curve Validation')
        plt.grid(True)

        plt.subplot(1,4,2)
        plt.plot(model_time_p_cc, model_p_cc, label='model output')
        plt.plot(exp_time_p_cc, exp_p_cc, label='data')
        plt.xlabel('Time (s)')
        plt.ylabel('Chamber Pressure (Pa)')
        plt.title('Chamber Pressure Validation')
        plt.grid(True)

        plt.subplot(1,4,3)
        plt.plot(model_time_p_ox_tank, model_p_ox_tank, label='model output')
        plt.plot(exp_time_p_ox_tank, exp_p_ox_tank, label='data')
        plt.xlabel('Time (s)')
        plt.ylabel('Ox Tank Pressure (Pa)')
        plt.title('Ox Tank Pressure Validation')
        plt.grid(True)

        plt.subplot(1,4,4)
        plt.plot(model_time_p_fuel_tank, model_p_fuel_tank, label='model output')
        plt.plot(exp_time_p_fuel_tank, exp_p_fuel_tank, label='experimental data')
        plt.xlabel('Time (s)')
        plt.ylabel('Fuel Tank Pressure (Pa)')
        plt.title('Fuel Tank Pressure Validation')
        plt.grid(True)

        plt.legend()
        plt.show()