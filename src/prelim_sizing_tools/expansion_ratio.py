from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel

import matplotlib

import numpy as np
import matplotlib.pyplot as plt

fuel_name = 'paraffin'
#C32H66 from RPA Paraffin Wax Composition
fuel_properties = f"""
fuel paraffin  C 32   H 66    wt%=100.00
h,KJ/Kgmol=-1860600     t(k)=298.15   rho,kg/m3={900}
"""

add_new_fuel(fuel_name, fuel_properties)

ox_name = 'N2O'

# Create CEA object for the first graph (Flame Temp/mw vs. Chamber Pressure)
C = CEA_Obj(
    oxName=ox_name, fuelName=fuel_name, pressure_units='Pa', isp_units='sec', cstar_units='m/s',
    temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='kJ/kg', density_units='kg/m^3',
    specific_heat_units='kJ/kg-K'
)


#mass flow rate input as constant
OF = 6
min_TW_ratio = 11
rocket_mass = 23 #kg


m_dot = 1.04 #kg/s
P_cc = 4e6 #Pa
#get cstar
c_star = C.get_Cstar(P_cc, OF)
#solve throat area using equation 3
A_throat = m_dot*c_star/P_cc
throat_diam = 39.3701*(2 *np.sqrt(A_throat/np.pi))

###setup array of atm pressures!

#assuming nasa common earth model
launch_pad_alt = 0 #(m) above sea level
apogee_height = 3048 #(m) from launch pad

height_arr = np.linspace(launch_pad_alt, (launch_pad_alt+apogee_height), 50)
pressure_arr = []
for h in height_arr:
    t = 15.04 - 0.00649*h #Celsius
    p = 1000 * 101.29*( (t+273.1)/288.08)**5.256 #Pa
    pressure_arr.append(p) #Pa

###pick ideal expansion altitude and solve expansion ratio
optimal_height = (2/3)*apogee_height #m above launch pad
P_exit = 1000 * 101.29*( ((15.04 - 0.00649*optimal_height)+273.1)/288.08)**5.256 #Pa

"""
first need to get ratio of specific heats (y), from testing noticed in range of expected chamber pressures
and expansion ratios y was constant to expratio and only changed with P_cc so we can solve here w/out losing accuracy
"""
fluid_prop = C.get_Chamber_MolWt_gamma(P_cc, OF, 5) #expratio guess 5
y = fluid_prop[1] # (-)

expratio = ( ( ((y+1)/2)**(1/(y-1)) ) * ( (P_exit/P_cc)**(1/y) ) * np.sqrt( ((y+1)/(y-1)) * ( (1- (P_exit/P_cc)**((y-1)/y) )) ) )**-1
print("expansion ratio: ", expratio)
A_exit = expratio * A_throat

thrust_arr = []
optimal_thrust_arr = []
#from an array of exit pressures (ideal expansion at some altitude)
for P_atm in pressure_arr:

    #solve equation 1 (need to call rocketcea to get ratio of specific heats)
    cf = np.sqrt( ((2*y**2)/(y-1)) * ( (2/(y+1))**((y+1)/(y-1)) ) * (1- (P_exit/P_cc)**((y-1)/y)) ) + (((P_exit-P_atm)*A_exit)/(P_cc*A_throat))
    optimal_cf = np.sqrt( ((2*y**2)/(y-1)) * ( (2/(y+1))**((y+1)/(y-1)) ) * (1- (P_atm/P_cc)**((y-1)/y)) )

    #solve equation 2 to get thrust
    thrust = cf * A_throat * P_cc
    optimal_thrust = optimal_cf * A_throat * P_cc

    #add thrust and convert exit pressure to altitude and add that to an array.
    thrust_arr.append(thrust)
    optimal_thrust_arr.append(optimal_thrust)

#plot vertical line for minimum thrust to weight ratio
min_start_thrust = (rocket_mass*9.81) * min_TW_ratio
plt.axvline(x=min_start_thrust, color='r', label=f'min starting thrust for T/W of {min_TW_ratio}')
#plot

plt.plot(thrust_arr, height_arr, label=f'optimal expansion at {optimal_height} (m)')
plt.plot(optimal_thrust_arr, height_arr, label=f'optimal expansion throughout burn')
plt.title(label=f'Altitude vs Thrust for throat_diam (in) {throat_diam}' )
plt.xlabel('Thrust (N)')
plt.ylabel('Altitude (m)')
plt.grid(True)
plt.legend()
plt.show()
