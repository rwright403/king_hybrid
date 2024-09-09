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
optimal_height = (2/3)*apogee_height #m above launch pad - #NOTE: this uses a rule of thumb, works rlly well for small change in alt sounding rocket flight
P_exit = 1000 * 101.29*( ((15.04 - 0.00649*optimal_height)+273.1)/288.08)**5.256 #Pa


OF_ratio = [4, 5, 6, 7, 8, 9]
chamber_pressures = [10, 20, 30, 40, 50]


###END USER INPUT

# Create CEA object for the first graph (Flame Temp/mw vs. Chamber Pressure)
ceaObj = CEA_Obj(
    oxName=ox_name, fuelName=fuel_name, pressure_units='Bar', isp_units='sec', cstar_units='m/s',
    temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='kJ/kg', density_units='kg/m^3',
    specific_heat_units='kJ/kg-K'
)


# Create subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 6))

# First subplot: Flame Temp/mw vs. Chamber Pressure
axs[0].set_title('Flame Temp/mw vs Chamber Pressure')
axs[0].set_xlabel('Chamber Pressure (Bar)')
axs[0].set_ylabel('Flame Temp/mw ((mol K)/kg)')

axs[1].set_title('Flame Temp vs Chamber Pressure')
axs[1].set_xlabel('Chamber Pressure (Bar)')
axs[1].set_ylabel('Flame Temp (K)')

for j in OF_ratio:
    flame_temp_mw_ratio_arr = []
    flame_temp_arr = []
    p_cc_arr = []

    for k in chamber_pressures:

        ###solve expansion ratio for each chamber pressure using y guess
        ###pick ideal expansion altitude and solve expansion ratio


        #first need to get ratio of specific heats (y), from testing noticed in range of expected chamber pressures
        #and expansion ratios y was constant to expratio and only changed with P_cc so we can solve here w/out losing accuracy

        #fluid_prop = CEA_Obj.get_Chamber_MolWt_gamma(k, j, 5.5) #expratio guess 5.5
        y = 1.23 #fluid_prop[1] # (-)

        expratio = ( ( ((y+1)/2)**(1/(y-1)) ) * ( (P_exit/k)**(1/y) ) * np.sqrt( ((y+1)/(y-1)) * ( (1- (P_exit/k)**((y-1)/y) )) ) )**-1


        i = ceaObj.get_IvacCstrTc_ChmMwGam(Pc=k, MR=j, eps=expratio)
        flame_temp_mw_ratio_arr.append(i[2] / i[3])
        flame_temp_arr.append(i[2])
        p_cc_arr.append(k)

    print(f"O/F Ratio: {j}, gamma = {i[4]}, gas const = {8134 / i[3]}, mw = {i[3]}")
    
    #graph OF vs flame temp / mw
    axs[0].plot(p_cc_arr, flame_temp_mw_ratio_arr, label=f'O/F={j}')
    #graph OF vs flame temp
    axs[1].plot(p_cc_arr, flame_temp_arr, label=f'O/F={j}')

axs[0].legend()
axs[0].grid(True)
axs[1].legend()
axs[1].grid(True)

# Second subplot: ISP vs. O/F Ratio at Different Chamber Pressures


axs[2].set_title('ISP vs O/F Ratio at Different Chamber Pressures')
axs[2].set_xlabel('O/F Ratio')
axs[2].set_ylabel('ISP (s)')

ispObj = CEA_Obj(
    oxName=ox_name, fuelName=fuel_name, pressure_units='Bar', isp_units='sec', cstar_units='m/s',
    temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='kJ/kg', density_units='kg/m^3',
    specific_heat_units='kJ/kg-K'
)

for Pc in chamber_pressures:


    of_arr = []
    isp_arr = []

    expratio = 5.5
    isp_prev = 0
    max_isp = 0
    of_store = 0
    i = 0.55

    while i < 16:
        isp = ispObj.get_Isp(Pc=Pc, MR=i, eps=expratio, frozen=0, frozenAtThroat=0)
        of_arr.append(i)
        isp_arr.append(isp)

        if isp > isp_prev:
            max_isp = isp
            of_store = i
        i = i + 0.5
        isp_prev = isp

    print(f"For Pc={Pc} bar, Max ISP: {max_isp}, O/F Ratio: {of_store}")
    #graph OF vs max theoretical ISP
    axs[2].plot(of_arr, isp_arr, label=f'Pc={Pc} bar [{bar_to_psi(Pc)} PSI]')

axs[2].legend()
axs[2].grid(True)

#show graphs, and TODO: print out stoichiometric O/F, O/F at highest flame temp/mw
plt.show()




selected_OF = 0
while(selected_OF == 0):
    selected_OF = float(input("Pick and Enter O/F ratio: "))

selected_Pcc = 0
while(selected_Pcc == 0):
    selected_Pcc = float(input("Pick and Enter P_cc (bar): "))




#selected_OF = 8
###NOTE: BUG: UNIT CONVERSION!!!!!
selected_Pcc *= 1e5



#step 2: solving the throat area

#use atmospheric model to get exit pressure at optimal height
#NOTE: DOCUMENT
P_pad = 1000 * 101.29*( ((15.04 - 0.00649*0)+273.1)/288.08)**5.256 #Pa #NOTE: sea level for getting base program
P_exit = 1000 * 101.29*( ((15.04 - 0.00649*optimal_height)+273.1)/288.08)**5.256 #Pa
P_apogee = 1000 * 101.29*( ((15.04 - 0.00649*apogee_height)+273.1)/288.08)**5.256 #Pa #NOTE: sea level for getting base program


#solve expansion ratio for optimal height
[print(P_exit, selected_Pcc)]
expratio = ( ( ((y+1)/2)**(1/(y-1)) ) * ( (P_exit/selected_Pcc)**(1/y) ) * np.sqrt( ((y+1)/(y-1)) * ( (1- (P_exit/selected_Pcc)**((y-1)/y) )) ) )**-1
print("expansion ratio: ", expratio)



#solve inital and final Cf, assuming rocket can always reach apogee (uses expansion ratio!!!!)
Cf_opt = np.sqrt( ((2*y**2)/(y-1)) * ( (2/(y+1))**((y+1)/(y-1)) ) * (1- (P_exit/selected_Pcc)**((y-1)/y)) ) 

#start calculate impulse with spreadsheet line of best fit eqn
It_est = 2.73*apogee_height + 4829

#use impulse estimation to graph a bunch of preliminary thrust curves that differ based on burn time and show the user
burn_time_arr = [2,3, 4, 5, 6, 7, 8, 9, 10]

for t in burn_time_arr:

    #solve throat area
    A_throat = It_est / (t*selected_Pcc*Cf_opt)

    #solve thrust
    F_thrust = selected_Pcc*A_throat*Cf_opt

    #plot curve
    
    time = [0, t, t]  # Example times (t0, t1, t2)
    thrust = [F_thrust, F_thrust, 0]  # Example thrust values (initial, final, 0)

    # Plot the thrust curve
    plt.plot(time, thrust, label=f't={t}, Throat Diam={39.3701*(np.sqrt(4*A_throat/np.pi))} (in)')

plt.title(label=f'Preliminary Thrust Curve Based on Estimated Total Impulse' )
plt.xlabel('Burn Time (s)')
plt.ylabel('Thrust (N)')
plt.grid(True)
plt.legend()
plt.show()


#allow the user to pick the burn time

selected_tburn = 0
while(selected_tburn == 0):
    selected_tburn = float(input("Pick and Enter Burn Time (s): "))



A_throat = It_est / (selected_tburn*selected_Pcc*Cf_opt)

#solve exit area

#solve exit velocity under those conditions

#solve mass flow rate #do we need this????

#solve Mprop




#ask user to update constants?





