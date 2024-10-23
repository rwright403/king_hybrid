CC_MODEL_MAP = {
    1: 'src.models.hybrid_cc_w_fuel_grain',
    2: 'src.models.adiabatic_lre_cc'
}

TANK_MODEL_MAP = {
    1: 'src.models.bens_ox_tank',
    2: 'src.models.adiabatic_ext_pressure_fed_cryo',
    3: 'src.models.adiabatic_pressurized_liquid_tank',
    4: 'src.models.piping_real_fuel_tank_data',
}


import matplotlib.pyplot as plt
import numpy as np
import csv
import importlib
import inspect

def get_variable_name(var):
    #Get the name of the variable as a string.
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    name = [var_name for var_name, var_val in callers_local_vars if var_val is var]
    return name[0] if name else "Unknown"

def to_csv(x_arr, y_arr, filename):
    combined_arr = list(zip(x_arr,y_arr))
    with open(r'./src/' + f'{filename}' + '.csv' , 'w', newline='') as file:
        #print(file)
        writer = csv.writer(file)
        writer.writerows(combined_arr)



def get_model(model_code, char):
    module_path = None
    if(char == "T"):
        module_path = TANK_MODEL_MAP.get(model_code, None)
    if(char == "C"):
        module_path = CC_MODEL_MAP.get(model_code, None)
 
    if module_path is None:
        raise ValueError(f"Model {model_code} is not defined.")
    
    module = importlib.import_module(module_path)
    return module.model  # or whichever class is appropriate for your tank model



def run_thrust_curve(inputs):

    ### SETUP FOR BOTH ENGINES:

    r1cc = None
    r1ox = None
    s1_fuel_tank = None

    time_arr = []
    m_dot_arr = []
    thrust_arr = []
    p_cc_arr = []
    p_ox_tank_arr = []
    p_fuel_tank_arr = []
    pressure_data = []
    
    #TODO: IS THIS THE BEST SPOT FOR THIS?? WORKS FOR NOW
    P_cc = inputs.P_atm

    ### Combustion Chamber Setup:
    cc_model_class = get_model(inputs.analysis_mode[0],'C')

    if inputs.analysis_mode[0] == 1:
        r1cc = cc_model_class(inputs.oxName, inputs.fuelName, inputs.CEA_fuel_str, inputs.m_fuel_i, 
                inputs.rho_fuel, inputs.a, inputs.n, inputs.L, inputs.A_port_i, 
                inputs.P_atm, inputs.A_throat, inputs.A_exit, inputs.TIMESTEP,)
        
    if inputs.analysis_mode[0] == 2:
        r1cc = cc_model_class(inputs.oxidizer_name, inputs.fuel_name, inputs.A_throat, inputs.A_exit, inputs.P_atm, inputs.TIMESTEP)


    ### Oxidizer Tank Setup:
    OxTank_model_class = get_model(inputs.analysis_mode[1], 'T')

    if inputs.analysis_mode[1] == 1:
        r1ox = OxTank_model_class(inputs.oxName, inputs.TIMESTEP, inputs.m_ox, inputs.Cd_1, inputs.A_inj_1,
                inputs.V_tank, inputs.P_tank, inputs.P_atm, inputs.all_error, inputs.inj_model)
    
    if inputs.analysis_mode[1] == 2:
        pass
        #TODO: how to handle two functions!!!!
        #pressurantTank = simpleAdiabaticExtPressurantTank(pressurant_name, P_prestank, 0.5, P_oxtank, 0.01, OUTLET_DIAM,TIMESTEP)
        #oxidizerTank = simpleAdiabaticPressurizedTank(oxidizer_name, pressurant_name, ID_PROPTANK, P_oxtank, 5, V_PROPTANK, TIMESTEP)
    
    if inputs.analysis_mode[1] == 3:
        r1ox = OxTank_model_class(inputs.pressurant_name, inputs.m_pressurant, inputs.fuel_name, inputs.m_fuel, inputs.P_fueltank, inputs.ID_PROPTANK, inputs.TIMESTEP)
        
    ### HYBRID THRUST CURVE
    if len(inputs.analysis_mode) == 2:
        pressure_data = [p_cc_arr, p_ox_tank_arr]
        #### HYBRID THRUST CURVE

        r1ox.inst(P_cc)
        #TODO: FIX with sim time? not sure its not working w OF now
        while (r1ox.t < 1):
            #print(r1cc.OF)
            r1cc.inst(r1ox.m_dot_ox)
            r1ox.inst(r1cc.P_cc)

            #add to arrays
            time_arr.append(r1ox.t)
            m_dot_arr.append(r1ox.m_dot_ox)
            thrust_arr.append(r1cc.instThrust)
            p_cc_arr.append(r1cc.P_cc) #NOTE: convert to atmospheric depending on data
            p_ox_tank_arr.append(r1ox.P_tank)
            #total_propellant = r1cc.total_propellant


            #print(r1ox.P_tank, r1cc.P_cc,r1ox.P_tank- r1cc.P_cc,  )

            #print(r1ox.m_ox, r1cc.m_fuel_t)

            #print(r1ox.t, r1cc.v_exit,r1cc.m_dot_cc_t,r1cc.R)

            #print("time: ", r1ox.t)

        #print("\n", "### mass balance ###")
        #print("total propellant calculated through sim: ", total_propellant+r1ox.m_ox+r1cc.m_fuel_t, "total starting/input propellant: ", inputs.m_ox+inputs.m_fuel_i, "difference (conservation of mass): ",total_propellant+r1ox.m_ox+r1cc.m_fuel_t -inputs.m_ox-inputs.m_fuel_i)


    ### for liquid setup fuel tank
    if len(inputs.analysis_mode) == 3:
        fuel_tank_model_class = get_model(inputs.analysis_mode[2],'T')

        pressure_data = [p_cc_arr, p_ox_tank_arr, p_fuel_tank_arr]

        if inputs.analysis_mode[2] == 1:
            print("model invalid for fuel tank (cannot use hybrid cc for liquid engine)")
        
        if inputs.analysis_mode[2] == 2:
            print("todo: implement")
            
        if inputs.analysis_mode[2] == 3:
            s1_fuel_tank = fuel_tank_model_class(inputs.pressurant_name, inputs.m_pressurant, inputs.fuel_name_1, inputs.m_fuel, inputs.P_fueltank, inputs.ID_PROPTANK, inputs.V_tank_2, inputs.Cd_2, inputs.A_inj_2, inputs.T_amb, inputs.TIMESTEP)
        
        if inputs.analysis_mode[2] == 4:
            s1_fuel_tank = fuel_tank_model_class(inputs.T_tank, inputs.Cd_spi, inputs.A_inj, inputs.fuel_name, inputs.fuel_tank_pressure_filepath)

        ### LIQUID THRUST CURVE
        
        time_arr = []
        m_dot_arr = []
        thrust_arr = []
        p_cc_arr = []
        p_ox_tank_arr = []
        p_fuel_tank_arr = []
        pressure_data = [p_cc_arr, p_ox_tank_arr, p_fuel_tank_arr]
        P_cc = inputs.P_atm

        ### Setup For Parsing ###

        ## Heat Transfer Parsing

        y_peak = 0
        cp_peak = 0
        P_cc_peak = 0
        C_star_peak = 0
        T_flame_peak = 0
        viscosity_peak = 0

        ## Injector Design Sheet Parsing

        # Oxidizer Inj
        smallest_ox_inj_pressure_drop = r1ox.P_tank
        m_dot_ox_min_dp = 0
        rho_ox_min_dp = 0
        kinematic_visc_ox_min_dp = 0
        y_ox_min_dp = 0
        t_ox_min_dp = 0
        p_ox_up_min_dp = 0
        p_ox_down_min_dp = 0
        #will also print: inputs.Cd_1

        # Fuel Inj
        smallest_fuel_inj_pressure_drop = s1_fuel_tank.P_tank
        m_dot_fuel_min_dp = 0
        rho_fuel_min_dp = 0
        kinematic_visc_fuel_min_dp = 0
        y_fuel_min_dp = 0
        t_fuel_min_dp = 0
        p_fuel_up_min_dp = 0
        p_fuel_down_min_dp = 0
        #will also print: inputs.Cd_2


        r1ox.inst(P_cc)
        s1_fuel_tank.inst(P_cc)
        

        m_fuel_burned = 0
        m_ox_burned = 0
        #TODO: FIX with sim time? not sure its not working w OF now
        while (r1ox.t < inputs.sim_time):
            #print(r1cc.OF)
            #BUG: cant handle 2 inputs to r1cc????
            #print("initial mass flow rates: ",r1ox.m_dot_ox, s1_fuel_tank.m_dot_fuel)
            #NOTE: fuel seems to be significantly underpredicting, and nan propagating through model
            r1cc.inst(r1ox.m_dot_ox, s1_fuel_tank.m_dot_fuel)
            #print("thrust curve: ",r1ox.t, r1cc.P_cc)
            r1ox.inst(r1cc.P_cc)
            s1_fuel_tank.inst(r1cc.P_cc)

            #print("looking at m_dot_f, m_dot_ox: ", s1_fuel_tank.m_dot_fuel, r1ox.m_dot_ox, r1cc.OF)
    
            m_fuel_burned += s1_fuel_tank.m_dot_fuel*r1ox.timestep
            m_ox_burned += r1ox.m_dot_ox*r1ox.timestep

            ### Parsing for detail design sheets ###NOTE: not super general if using other models, if it breaks in future read this and comment v out
            if (r1cc.instThrust > r1cc.prev_thrust):

                
                arr1 = r1cc.C.get_Throat_MolWt_gamma(r1cc.P_cc, r1cc.OF, r1cc.expratio, frozen=0)
                arr2 = r1cc.C.get_Temperatures(r1cc.P_cc, r1cc.OF, r1cc.expratio, frozen=0, frozenAtThroat=0)
                arr3 = r1cc.C.get_Throat_Transport(r1cc.P_cc, r1cc.OF, r1cc.expratio, frozen=0)
                
                y_peak = arr1[1]
                cp_peak = arr3[0]
                P_cc_peak = r1cc.P_cc*(2/(y_peak+1))**(y_peak/(y_peak-1))
                C_star_peak = r1cc.C.get_Cstar(r1cc.P_cc, r1cc.OF)
                T_flame_peak = arr2[1]
                viscosity_peak = arr3[1]
                R_peak = r1cc.R

            inst_ox_inj_pressure_drop = r1ox.P_tank - r1cc.P_cc
            if (inst_ox_inj_pressure_drop < smallest_ox_inj_pressure_drop):

                smallest_ox_inj_pressure_drop = inst_ox_inj_pressure_drop

                p_ox_up_min_dp = r1ox.P_tank
                p_ox_down_min_dp = r1cc.P_cc
                
                m_dot_ox_min_dp = r1ox.m_dot_ox
                rho_ox_min_dp = r1ox.rho_exit
                #kinematic_visc_ox_min_dp = r1ox.kinematic_visc_ox
                #y_ox_min_dp = r1ox.y_ox
                
                t_ox_min_dp = r1ox.t



            inst_fuel_inj_pressure_drop = s1_fuel_tank.P_tank - r1cc.P_cc
            if (inst_fuel_inj_pressure_drop < smallest_fuel_inj_pressure_drop):

                smallest_fuel_inj_pressure_drop = inst_fuel_inj_pressure_drop

                p_fuel_up_min_dp = s1_fuel_tank.P_tank
                p_fuel_down_min_dp = r1cc.P_cc

                m_dot_fuel_min_dp = s1_fuel_tank.m_dot_fuel
                rho_fuel_min_dp = s1_fuel_tank.rho_prop
                kinematic_visc_fuel_min_dp = s1_fuel_tank.kinematic_visc_fuel
                y_fuel_min_dp = s1_fuel_tank.y_fuel

                t_fuel_min_dp = r1ox.t



            #RECORD DATA
            time_arr.append(r1ox.t)
            m_dot_arr.append(r1cc.m_dot_cc_t)
            thrust_arr.append(r1cc.instThrust)
            p_cc_arr.append(r1cc.P_cc)
            p_ox_tank_arr.append(r1ox.P_tank)
            p_fuel_tank_arr.append(s1_fuel_tank.P_tank)

            #print(s1_fuel_tank.m_fuel)

            #print(r1cc.instThrust,s1_fuel_tank.P_tank)
            
        #print("total propellant burned in simshould match report: ", m_fuel_burned, m_ox_burned)
            
    total_impulse = np.trapz(thrust_arr, time_arr)
    print(f"Total Engine Impulse: {total_impulse:.6f} (N s)")

    ###WRITE CSV FOR FLIGHT SIM, VALIDATION AND OTHER EXTERNAL ANALYSIS
    to_csv(time_arr,m_dot_arr, "m_dot_cc_t")
    to_csv(time_arr,thrust_arr, "thrust")
    to_csv(time_arr,p_cc_arr, "p_cc")
    to_csv(time_arr,p_ox_tank_arr, "p_ox_tank")
    to_csv(time_arr,p_fuel_tank_arr, "p_fuel_tank")


    ###PLOTS
    if inputs.thrust_curve_graphs == True:

        plt.subplot(1,2,1)
        plt.plot(time_arr,thrust_arr)
        plt.xlabel('Time (s)')
        plt.ylabel('Thrust (N)')
        plt.title('Thrust Curve')
        plt.grid(True)

        plt.subplot(1,2,2)
        for i in pressure_data:
            plt.plot(time_arr,i,label=get_variable_name(i))
        plt.xlabel('Time (s)')
        plt.ylabel('Pressure (Pa)')
        plt.title('System Pressures Over Time')
        plt.grid(True)
        plt.legend()


        print(f"\nThroat Properties at Peak Thrust for Heat Transfer\n------------\nRatio of specific heats: {y_peak} (-)\nSpec. Heat Const. Pres. {cp_peak} (J/(kg K))\nThroat Pressure {P_cc_peak} (Pa)\nCharacteristic Velocity {C_star_peak} (m/s)\nThroat Flame Temp {T_flame_peak} (K)\nViscosity {viscosity_peak} (Pa s)\nGas Constant {R_peak} (J/(kg K))")

        print(f"\nMinimum Pressure Drop Fuel Inj Properties for Sizing\n------------\nTotal Fuel Mass Flow rate of all elements: {m_dot_fuel_min_dp} (kg/s)\nUpstream Pressure at inst: {p_fuel_up_min_dp} (Pa)\nDownstream Pressure at inst: {p_fuel_down_min_dp} (Pa)\nFuel Density at Orifice Outlet {rho_fuel_min_dp} (kg/m^3)\nFuel Kinematic Viscosity {kinematic_visc_fuel_min_dp} (Pa s)\nFuel Ratio of specific heats: {y_fuel_min_dp} (-)\nFuel Orifice Discharge Coeff: {inputs.Cd_2} (-)\nAt t = {t_fuel_min_dp} (s)")
        print(f"\nMinimum Pressure Drop Ox Inj Properties for Sizing\n------------\nTotal Ox Mass Flow rate of all elements: {m_dot_ox_min_dp} (kg/s)\nUpstream Pressure at inst: {p_ox_up_min_dp} (Pa)\nDownstream Pressure at inst: {p_ox_down_min_dp} (Pa)\nOx Density at Orifice Outlet {rho_ox_min_dp} (kg/m^3)\nOx Kinematic Viscosity {kinematic_visc_ox_min_dp} (Pa s)\nOx Ratio of specific heats: {y_ox_min_dp} (-)\nOx Orifice Discharge Coeff: {inputs.Cd_1} (-)\nAt t = {t_ox_min_dp} (s)")

        plt.show()
