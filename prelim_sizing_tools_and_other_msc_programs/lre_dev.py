from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel

import matplotlib
matplotlib.use('Qt5Agg')  # Use the TkAgg backend

import numpy as np
import matplotlib.pyplot as plt


def bar_to_psi(x):
    return x * 14.5038

###START USER INPUT
fuel_name = 'Propane'

#fuel_name = 'paraffin'
##C32H66 from RPA Paraffin Wax Composition
#fuel_properties = f"""
#fuel paraffin  C 32   H 66    wt%=100.00
#h,KJ/Kgmol=-1860600     t(k)=298.15   rho,kg/m3={900}
#"""

#add_new_fuel(fuel_name, fuel_properties)


ox_name = 'N2O'

apogee_height = 3048 #m
optimal_height = (2/3)*apogee_height #m above launch pad
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

    #NOTE: bad var name!
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
    axs[2].plot(of_arr, isp_arr, label=f'Pc={Pc} bar  [{bar_to_psi(Pc)} PSI]')

axs[2].legend()
axs[2].grid(True)

#show graphs, and TODO: print out stoichiometric O/F, O/F at highest flame temp/mw
plt.show()

### step 2 ###








"""
#note expratio guess
#expratio = 5
"""

selected_OF = 0
while(selected_OF == 0):
    selected_OF = float(input("Pick and Enter O/F ratio: "))

selected_Pcc = 0
while(selected_Pcc == 0):
    selected_Pcc = float(input("Pick and Enter P_cc (bar): "))
"""
"""
selected_OF = 8
selected_Pcc = 40

i = ceaObj.get_IvacCstrTc_ChmMwGam(selected_Pcc, selected_OF, expratio)
Isp = i[0]
T_flame = i[2] 
Mw = i[3]
y = i[4]


### plot CF vs P1/P2 ###
P_atm = 1 #bar
expratio = [4,5,10,80,100,200]
P_exit = [0.7,0.8,0.9,1,1.1,1.2]

#note: we have selected P_cc so y will be constant!!!!!

for l in expratio:
    test = ceaObj.get_IvacCstrTc_ChmMwGam(selected_Pcc, selected_OF, l)
    print("next expratio")

    P_ratio_arr = []
    Cf_arr = []


    for m in P_exit:

        Cf = np.sqrt( ((2*y**2)/(y-1))*(2/(y+1))**((y+1)/(y-1))*(1-(m/selected_Pcc)**((y-1)/y)) ) + ((m-P_atm)/selected_Pcc)*l
        P_ratio_arr.append(selected_Pcc/m)
        Cf_arr.append(Cf)

        print(((m-P_atm)/selected_Pcc)*l)
    
    plt.plot(P_ratio_arr, Cf_arr, label=f'expratio={l}')


#now we know that design will deviate from selected P_cc through operation and iteration so want to look at impact of P_cc change 
#abstract P_cc to expratio array and use gamma for each line!!!!!

#for a given expansion ratio
expratio = 1



plt.title(label=f'Cf vs P_cc/P_exit')
plt.xlabel('P_cc/P_exit')
plt.ylabel('Cf')
plt.legend()
plt.grid(True)
plt.show()

"""
input_of = 0
while(input_of ==0):
    input_of = input("Pick and Enter O/F ratio: ")

#now graph chamber pressure vs m dot for a given A_throat once selected O/F ratio
    


P_cc_min = chamber_pressures[0] #Bar
P_cc_max = chamber_pressures[-1] #Bar
P_cc_step = 2 #Bar
P_cc_current = P_cc_min

throat_diam_min = 0.75 #in
throat_diam_max = 1.75 #in
throat_diam_step = 0.125 #in
throat_diam_current = throat_diam_min

while(throat_diam_current<=throat_diam_max):
    mdot_arr = []
    p_cc_arr.clear()
    #solve throat area
    A_throat = (np.pi/4)*((0.0254*throat_diam_current)**2)

    while(P_cc_current<=P_cc_max):
                #call cea object - solve cstar
        cstar = ceaObj.get_Cstar(P_cc_current, input_of)
        print(cstar)

        #solve mdot / A throat
        j = P_cc_current/cstar
        #add to arrays
        mdot_arr.append(cstar)
        p_cc_arr.append(P_cc_current)

        P_cc_current+=P_cc_step
        #print(cstar)
    plt.plot(p_cc_arr, mdot_arr, label=f'throat area (m^2) {A_throat:.6f}')
    P_cc_current=P_cc_min
    throat_diam_current+=throat_diam_step


plt.title(label=f'mdot vs pressure for '+ ox_name + ' ' + fuel_name + f' at O/F Ratio {input_of}')
plt.xlabel('Chamber Pressure (Bar)')
plt.ylabel('mdot (kg/s)')
plt.legend()
plt.grid(True)
plt.show()
"""