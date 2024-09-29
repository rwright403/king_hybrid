#estimate vapor pressure wrt conditions
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

from rocketprops.rocket_prop import get_prop
from rocketprops.line_supt import calc_line_vel_dp
from rocketprops.valve_supt import cv_valve_dp


# Define temperature range (input in celsius, program converts to kelvin)
T_min = 16  # Minimum operating temperature
T_max = 27  # Maximum operating temperature
T_crit = 36 # Critical temperature of nitrous
num_points = 100  # Number of points for the plot

#INPUT FOR FEED SYSTEM:
m_dot = 1.14 #kg/s
#TODO: CHECK IF THIS IS CORRECT***** GO THROUGH PRELIM DESIGN SS

id_arr = [0.25, 0.375, 0.5] #id of pipe you are considering (in))

first_len = 12 #in
first_k_fac_arr = np.array([0.4]) #input all K factors for pipe between tank and valve
first_rough = 4e-7 #roughness (in)

cv = 20 #imperial cv for valve

second_len = 12 #in
second_k_fac_arr = np.array([0]) #input all K factors for pipe between valve and injector
second_rough = 4e-7 #roughness (in)


# Generate an array of temperatures
T_min = T_min + 273.15
T_max = T_max + 273.15
T_crit = T_crit + 273.15
temperatures = np.linspace(T_min, T_crit, num_points)
temperatures_celsius = temperatures - 273.15  # Convert from Kelvin to Celsius
vapor_pressures = np.zeros_like(temperatures)

# Plotting
plt.figure(figsize=(10, 6))

# Calculate vapor pressure for each temperature using CoolProp
for i, T in enumerate(temperatures):
  P = CP.PropsSI('P', 'T', T, 'Q', 0, 'N2O')
  vapor_pressures[i] = P / 100000.0  # Convert from Pa to bar

# Plot vapor pressure for the current quality
plt.plot(temperatures_celsius, vapor_pressures, color='k', label='N2O Saturated Liquid Line')

# THESE TEMPERATURES ARE MAX AND MIN VALUES FOR HOTFIRING
plt.axvline(x=T_min-273.15, color='b', linestyle='--', label=f'{T_min-273.15} °C')
plt.axvline(x=T_max-273.15, color='c', linestyle='--', label=f'{T_max-273.15} °C')
plt.axhline(y=72.5, color='r', linestyle='--', label="N2O Critical Pressure")

plt.xlabel('Temperature (°C)')
plt.ylabel('Vapor Pressure (bar)')
plt.title('Vapor Pressure of Nitrous Oxide (N2O) vs. Temperature')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

#adjust this depending on what pressure of tank you want for feed system analysis
P_tank = np.interp(T_min-273.15, temperatures_celsius, vapor_pressures)

#rocketprops setup:
pObj = get_prop('N2O') #setup propellant object
T_op_R = 1.8*T_min #convert temp to Rankine
P_tank = P_tank * 14.5038 #convert bar to psia
m_dot = m_dot * 2.205 #convert kg/s to lb/s

first_k_fac = np.sum(first_k_fac_arr)
second_k_fac = np.sum(second_k_fac_arr)

#graph setup
x_arr = [0,first_len, first_len, first_len+second_len] #in

p_inj_arr = []

for id in id_arr:
  p_arr = [P_tank/14.5038] #bar

  #first line vel dp
  first_velo, first_deltaP = calc_line_vel_dp( pObj, T_op_R, P_tank, m_dot, id, first_rough, first_k_fac, first_len)
  p_valve_inlet = P_tank - first_deltaP
  p_arr.append(p_valve_inlet/14.5038) #bar

  #input pressure into cv valve dp
  valve_deltaP = cv_valve_dp( pObj, cv, m_dot, T_op_R, p_valve_inlet)
  p_second_inlet = p_valve_inlet - valve_deltaP
  p_arr.append(p_second_inlet/14.5038) #bar

  #second line ve dp
  second_velo, second_deltaP = calc_line_vel_dp( pObj, T_op_R, p_second_inlet, m_dot, id, second_rough, second_k_fac, second_len)
  p_inj_inlet = p_second_inlet - second_deltaP
  p_arr.append(p_inj_inlet/14.5038) #bar
  p_inj_arr.append(p_inj_inlet/14.5038) #bar

  #graph:
  plt.plot(x_arr, p_arr, label={id})

  print(first_velo, second_velo)

plt.xlabel('distance along feed (in)')
plt.ylabel('Pressure (bar)')
plt.title(f'Pressure drop (bar) along feed line (in) at {T_min-273.15} °C')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

for i in p_inj_arr:
  print(i)