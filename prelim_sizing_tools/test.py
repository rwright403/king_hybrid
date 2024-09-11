from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel

import matplotlib
matplotlib.use('Qt5Agg')  # Use the TkAgg backend

import numpy as np
import matplotlib.pyplot as plt


def bar_to_psi(x):
    return x * 14.5038

###START USER INPUT

fuel_name = 'paraffin'
#C32H66 from RPA Paraffin Wax Composition
fuel_properties = f"""
fuel paraffin  C 32   H 66    wt%=100.00
h,KJ/Kgmol=-1860600     t(k)=298.15   rho,kg/m3={900}
"""

add_new_fuel(fuel_name, fuel_properties)


ox_name = 'N2O'

apogee_height = 3048 #m
optimal_height = (2/3)*apogee_height #m above launch pad


###END USER INPUT

# Create CEA object for the first graph (Flame Temp/mw vs. Chamber Pressure)
ceaObj = CEA_Obj(
    oxName=ox_name, fuelName=fuel_name, pressure_units='Pa', isp_units='sec', cstar_units='m/s',
    temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='kJ/kg', density_units='kg/m^3',
    specific_heat_units='kJ/kg-K'
)


#guess impulse total
I_t = 10000 #N s

#troposphere model valid from 0m to 10,000m - https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html 
def troposphere_model(h):
    t = 15.04 - 0.00649*h #Celsius
    p = 1000 * 101.29*( (t+273.1)/288.08)**5.256 #Pa
    return p

#inputs

expratio = 5 #init guess
selected_OF = 8
selected_Pcc = 40e5
P_atm =1e5

i = ceaObj.get_IvacCstrTc_ChmMwGam(selected_Pcc, selected_OF, expratio)
Isp = i[0]
T_flame = i[2] 
Mw = i[3]
y = i[4]

apogee_height


P_pad = troposphere_model(0)
P_opt = troposphere_model(optimal_height)
P_apo = troposphere_model(apogee_height)


Cf_apo = np.sqrt( ((2*y**2)/(y-1))*(2/(y+1))**((y+1)/(y-1))*(1-(P_apo/selected_Pcc)**((y-1)/y)) ) + ((P_apo-P_atm)/selected_Pcc)*expratio
Cf_pad = np.sqrt( ((2*y**2)/(y-1))*(2/(y+1))**((y+1)/(y-1))*(1-(P_pad/selected_Pcc)**((y-1)/y)) ) + ((P_pad-P_atm)/selected_Pcc)*expratio

Cf_opt = np.sqrt( ((2*y**2)/(y-1))*(2/(y+1))**((y+1)/(y-1))*(1-(P_opt/selected_Pcc)**((y-1)/y)) ) + ((P_opt-P_atm)/selected_Pcc)*expratio

burn_times = np.linspace(1, 10, 10)

for t in burn_times:
    A_throat = 2*I_t / (t*selected_Pcc*(Cf_apo+Cf_pad))
    print(t, A_throat, 39.3701*np.sqrt(4*A_throat/np.pi), selected_Pcc*A_throat*Cf_opt)