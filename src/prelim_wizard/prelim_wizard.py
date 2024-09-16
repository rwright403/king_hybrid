from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel
import CoolProp.CoolProp as CP #I love coolprop! ~ units: http://www.coolprop.org/v4/apidoc/CoolProp.html

import matplotlib
matplotlib.use('Qt5Agg')  # Use the TkAgg backend

import numpy as np
import matplotlib.pyplot as plt


def bar_to_psi(x):
    return x * 14.5038


def magic(inputs):

    ### Setup
    apogee_height = 3048 #m
    optimal_height = (2/3)*apogee_height #m above launch pad - #NOTE: this uses a rule of thumb, works rlly well for small change in alt sounding rocket flight
    P_exit = 1000 * 101.29*( ((15.04 - 0.00649*optimal_height)+273.1)/288.08)**5.256 #Pa


    OF_ratio = [4, 5, 6, 7, 8, 9]
    chamber_pressures = [10, 20, 30, 40, 50]


    ### Step 1: setup Pcc, O/F, and combustion property plots to select P_cc and O/F

    # Create CEA object for the first graph (Flame Temp/mw vs. Chamber Pressure)
    ceaObj = CEA_Obj(
        oxName=inputs.oxidizer_name, fuelName=inputs.fuel_name, pressure_units='Bar', isp_units='sec', cstar_units='m/s',
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

    C = CEA_Obj(
        oxName=inputs.oxidizer_name, fuelName=inputs.fuel_name, pressure_units='Pa', isp_units='sec', cstar_units='m/s',
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
            isp = C.get_Isp(Pc=(Pc)*1e5, MR=i, eps=expratio, frozen=0, frozenAtThroat=0)
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

    print("Pick O/F and P_cc, close graph and enter when prompted:")
    plt.show()

    selected_OF = 0
    while(selected_OF == 0):
        selected_OF = float(input("Enter O/F ratio: "))

    selected_Pcc = 0
    while(selected_Pcc == 0):
        selected_Pcc = float(input("Enter P_cc (bar): "))




    ###NOTE: UNIT CONVERSION FROM BAR TO Pa!!!!!
    selected_Pcc *= 1e5

    ### Step 2: solving the throat area

    #use atmospheric model to get exit pressure at optimal height
    #NOTE: DOCUMENT

    P_pad = 1000 * 101.29*( ((15.04 - 0.00649*0)+273.1)/288.08)**5.256 #Pa #NOTE: sea level for getting base program
    P_exit = 1000 * 101.29*( ((15.04 - 0.00649*optimal_height)+273.1)/288.08)**5.256 #Pa
    P_apogee = 1000 * 101.29*( ((15.04 - 0.00649*apogee_height)+273.1)/288.08)**5.256 #Pa #NOTE: sea level for getting base program


    #solve expansion ratio for optimal height
    #print(P_exit, selected_Pcc)
    expratio = ( ( ((y+1)/2)**(1/(y-1)) ) * ( (P_exit/selected_Pcc)**(1/y) ) * np.sqrt( ((y+1)/(y-1)) * ( (1- (P_exit/selected_Pcc)**((y-1)/y) )) ) )**-1
    #print("expansion ratio: ", expratio)


    #solve inital and final Cf, assuming rocket can always reach apogee (uses expansion ratio!!!!)
    Cf_opt = np.sqrt( ((2*y**2)/(y-1)) * ( (2/(y+1))**((y+1)/(y-1)) ) * (1- (P_exit/selected_Pcc)**((y-1)/y)) ) 

    #start calculate impulse and rocket dry mass with spreadsheet line of best fit eqn
    It_est = 2.73*apogee_height + 4829
    rocket_dry_mass = (1.03e-3)*It_est + 21

    #use impulse estimation to graph a bunch of preliminary thrust curves that differ based on burn time and show the user
    burn_time_arr = [2, 3, 4, 5, 6, 7, 8, 9, 10]

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

    print("The following solved by assuming ideal expansion (calculated at 2/3rds altitude) occurs throughout burn.\nThis is a good assumption for low altitude suborbital rockets")

    min_start_thrust = (rocket_dry_mass*9.81) * inputs.min_TW_ratio
    print(f"For reference: to achieve the estimated min_start_thrust of {inputs.min_TW_ratio}, need a starting thrust of {min_start_thrust:.2f} (N).")

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


    ### Step 3: Evaluate Nozzle Performance
    A_throat = It_est / (selected_tburn*selected_Pcc*Cf_opt)
    throat_diam = 39.3701*(2 *np.sqrt(A_throat/np.pi))

    #solve exit area
    expratio = ( ( ((y+1)/2)**(1/(y-1)) ) * ( (P_exit/selected_Pcc)**(1/y) ) * np.sqrt( ((y+1)/(y-1)) * ( (1- (P_exit/selected_Pcc)**((y-1)/y) )) ) )**-1
    print("expansion ratio: ", expratio)
    A_exit = expratio * A_throat


    height_arr = np.linspace(inputs.elevation, (inputs.elevation+apogee_height), 50)
    pressure_arr = []
    for h in height_arr:
        t = 15.04 - 0.00649*h #Celsius
        p = 1000 * 101.29*( (t+273.1)/288.08)**5.256 #Pa
        pressure_arr.append(p) #Pa

    thrust_arr = []
    optimal_thrust_arr = []

    for P_atm in pressure_arr:

        #solve equation 1 (need to call rocketcea to get ratio of specific heats)
        cf = np.sqrt( ((2*y**2)/(y-1)) * ( (2/(y+1))**((y+1)/(y-1)) ) * (1- (P_exit/selected_Pcc)**((y-1)/y)) ) + (((P_exit-P_atm)*A_exit)/(selected_Pcc*A_throat))
        optimal_cf = np.sqrt( ((2*y**2)/(y-1)) * ( (2/(y+1))**((y+1)/(y-1)) ) * (1- (P_atm/selected_Pcc)**((y-1)/y)) )

        #solve equation 2 to get thrust
        thrust = cf * A_throat * selected_Pcc
        optimal_thrust = optimal_cf * A_throat * selected_Pcc

        #add thrust and convert exit pressure to altitude and add that to an array.
        thrust_arr.append(thrust)
        optimal_thrust_arr.append(optimal_thrust)

    #plot vertical line for minimum thrust to weight ratio
    min_start_thrust = (rocket_dry_mass*9.81) * inputs.min_TW_ratio
    plt.axvline(x=min_start_thrust, color='r', label=f'min starting thrust for T/W of {inputs.min_TW_ratio}')
    #plot

    print("Summary of Nozzle Performance:")
    plt.plot(thrust_arr, height_arr, label=f'optimal expansion at {optimal_height} (m)')
    plt.plot(optimal_thrust_arr, height_arr, label=f'optimal expansion throughout burn')
    plt.title(label=f'Altitude vs Thrust for throat_diam (in) {throat_diam}' )
    plt.xlabel('Thrust (N)')
    plt.ylabel('Altitude (m)')
    plt.grid(True)
    plt.legend()
    plt.show()


    #solve exit velocity under those conditions
    fluid_prop = C.get_Chamber_MolWt_gamma(selected_Pcc, selected_OF, expratio)
    R = 8314 / fluid_prop[0] # J/(kg K)
    y = fluid_prop[1] # (-)

    temperatures = C.get_Temperatures(selected_Pcc, selected_OF, expratio, 0, 1)
    T_cc = temperatures[0]

    v_exit = np.sqrt(((2 * y) / (y - 1)) * R * T_cc * (1 - (P_exit / selected_Pcc)**((y - 1) / y)))

    #TODO: solve mass flow rate #do we need this???? --> steady state performance at an instant is nice
    m_dot_cc = thrust / v_exit

    m_dot_fuel = m_dot_cc/(selected_OF +1)
    m_dot_ox = m_dot_cc - m_dot_fuel

    #pick tank pressures:
    selected_P_ox_tank = 0
    while(selected_P_ox_tank == 0):
        selected_P_ox_tank = float(input("Pick and Enter Oxidizer Tank Pressure (Bar): "))*1e5

    selected_P_fuel_tank = 0
    while(selected_P_fuel_tank == 0):
        selected_P_fuel_tank = float(input("Pick and Enter Fuel Tank Pressure (Bar): "))*1e5


    #estimate injector area!!!!
    print(f"estimating reqiured injector areas for fuel and oxidizer assuming a Cd of {inputs.Cd_est} and SPI model")
    
    #estimate propellant densities with coolprop
    rho_fuel = CP.PropsSI('D', 'P', selected_P_fuel_tank, 'T', 275, inputs.fuel_name)
    rho_oxidizer = CP.PropsSI('D', 'P', selected_P_fuel_tank, 'T', 275, inputs.oxidizer_name)
    
    A_ox_inj = m_dot_ox / (inputs.Cd_est*np.sqrt(2*rho_oxidizer*(selected_P_ox_tank-selected_Pcc)))

    A_fuel_inj = m_dot_fuel / (inputs.Cd_est*np.sqrt(2*rho_fuel*(selected_P_fuel_tank-selected_Pcc)))
    
    ### Finally, print a summary of the results

    print(f"\n\n\n------------\nPreliminary Design Summary:\n------------\nCombustion Chamber:\nP_cc: {selected_Pcc} (Pa)\nO/F Ratio: {selected_OF}")
    print(f"Flame Temperature: {T_cc}\nRatio of Specific Heats: {y}\nReactant Gas Const. {R} (J/(kg K))")
    print(f"Total Impulse: {It_est} (N s)\nAverage Thrust {thrust} (N)\nMass Flow Rate {m_dot_cc} (kg/s)\nExit Velocity: {v_exit} (m/s)\nBurn Time (s)")
    print(f"Expansion Ratio: {expratio}\nThroat Diameter {throat_diam} (in)")
    print(f"Feed System Pressures:\nOxidizer Tank Pressure: {selected_P_ox_tank} (Pa)\nFuel Tank Pressure: {selected_P_fuel_tank} (Pa)")
    print(f"Preliminary Injector Solved by Assuming Fuel and Oxidizer Orifices Follow SPI Model and a Cd guess of {inputs.Cd_est}")
    print(f"Oxidizer Injector Discharge Area: {A_ox_inj} (m^2), Fuel Injector Discharge Area {A_fuel_inj} (m^2)")
    print("------------\n")
    confirm = 0
    while(confirm == 0):
        confirm = float(input("Copy Results into input file? (Confirm -->'Y' /Any Other Char to Exit)"))

        if((confirm != 'Y') or (confirm != 'y')):
            exit()

        else:
            #TODO: Implement
            print("NOT IMPLEMENTED YET L MY BAD GUYS JUST COPY PASTE FOR NOW")

            """ 2 --> adiabatic_lre_cc"""
            #

            oxidizer_name = oxidizer_name
            fuel_name = fuel_name
            A_throat = 0.00102028641 #m^2
            A_exit = 0.00527334324 #m^2
            P_atm = P_atm
            TIMESTEP = timestep

            """ 1 --> bens_ox_tank"""
            #

            oxName = oxidizer_name
            timestep = timestep 
            m_ox = 4.48 #kg 
            #NOTE: GUESSING Cd
            C_inj_1 =  0.35 * 0.00007471705 #1* 0.00001735222#(num_orifices * Cd * orifice_diam) Note: guessing Cd of 0.6, NOTE: when it doesnt work this is why :)
            V_tank = 6.4e-3 # - from report: "5.8L of nos in a 6.4L tank"
            P_tank = 5.171e6 #Pa
            P_atm = P_atm 
            all_error = 0.01

            """simpleAdiabaticPressurizedTank"""
            #

            pressurant_name = pressurant_name
            m_pressurant  = 0.12 #NOTE: estimated for now based on volume they gave in report, should i change inputs to this model?
            fuel_name = fuel_name #NOTE: This might not work, assuming 100% when they used 95% as well
            m_fuel = 1.12 #kg 
            P_fueltank = 4.82633e6 #Pa
            ID_PROPTANK = 0.0254*5 #m 
            V_tank_2 = 2.16e-3 #m^3
            C_inj_2 = 0.6*0.0000136284 #m^2
            T_amb = T_amb
            TIMESTEP = timestep


    #update variables in input file!!!!









