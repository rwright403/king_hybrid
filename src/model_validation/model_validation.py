import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from src import constants

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

###READ FROM FILES!!!!
model_time_thrust, model_thrust = read_csv(constants.model_thrust_file_path)
model_time_p_cc, model_p_cc = read_csv(constants.model_p_cc_file_path)
model_time_p_tank, model_p_tank = read_csv(constants.model_p_tank_file_path)

exp_time_thrust, exp_thrust = read_csv(constants.exp_thrust_file_path)
exp_time_p_cc, exp_p_cc = read_csv(constants.exp_p_cc_file_path)
exp_time_p_tank, exp_p_tank = read_csv(constants.exp_p_tank_file_path)

#return summary of rocket!!!!
print("\n")
print("Rocket: ", constants.name)
print("\n","### injector model: ###")
if(constants.inj_model == 1):
    print("SPI")
elif(constants.inj_model == 2):
    print("HEM")
elif(constants.inj_model == 3):
    print("DYER")
print("\n")

###PLOT!!!
plt.subplot(1,3,1)
plt.plot(exp_time_thrust, exp_thrust,label='data') 
plt.plot(model_time_thrust, model_thrust, label='model output')
plt.xlabel('Time (s)')
plt.ylabel('Thrust (N)')
plt.title('Thrust Curve Validation')
plt.grid(True)

plt.subplot(1,3,2)
plt.plot(exp_time_p_cc, exp_p_cc, label='data')
plt.plot(model_time_p_cc, model_p_cc, label='model output')
plt.xlabel('Time (s)')
plt.ylabel('Chamber Pressure (Pa)')
plt.title('Chamber Pressure Validation')
plt.grid(True)

plt.subplot(1,3,3)
plt.plot(exp_time_p_tank, exp_p_tank, label='data')
plt.plot(model_time_p_tank, model_p_tank, label='model output')
plt.xlabel('Time (s)')
plt.ylabel('Tank Pressure (Pa)')
plt.title('Tank Pressure Validation')
plt.grid(True)

plt.legend()
plt.show()