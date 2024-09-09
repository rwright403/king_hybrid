CC_MODEL_MAP = {
    1: 'src.models.hybrid_cc_w_fuel_grain',
    2: 'adiabatic_lre_cc'
}

TANK_MODEL_MAP = {
    1: 'src.models.bens_ox_tank',
    2: 'src.models.adiabatic_ext_pressure_fed_cryo',
    3: 'src.tank_models.adiabatic_pressurized'
}


import matplotlib.pyplot as plt
import numpy as np
import csv
import importlib

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


    print(str,module_path)    
    if module_path is None:
        raise ValueError(f"Tank model {model_code} is not defined.")
    
    module = importlib.import_module(module_path)
    return module.model  # or whichever class is appropriate for your tank model



def run_thrust_curve(inputs):

    ### HYBRID SETUP
    if len(inputs.analysis_mode) == 2:

        #input desired hybrid cc w fuel grain model and setup:
        cc_model_class = get_model(inputs.analysis_mode[0],'C')


        r1cc = cc_model_class(inputs.oxName, inputs.fuelName, inputs.CEA_fuel_str, inputs.m_fuel_i, 
                inputs.rho_fuel, inputs.a, inputs.n, inputs.L, inputs.A_port_i, 
                inputs.P_atm, inputs.A_throat, inputs.A_exit, inputs.timestep,)

        #input ox tank model and setup:
        OxTank_model_class = get_model(inputs.analysis_mode[1], 'T')

        r1ox = OxTank_model_class(inputs.oxName, inputs.timestep, inputs.m_ox, inputs.C_inj,
                inputs.V_tank, inputs.P_tank, inputs.P_atm, inputs.all_error, inputs.inj_model)
        
        #### HYBRID THRUST CURVE
        time_arr = []
        m_dot_arr = []
        thrust_arr = []
        p_cc_arr = []
        p_tank_arr = []
        P_cc = inputs.P_atm

        r1ox.inst(P_cc)
        while (r1cc.OF != 0):
            #print(r1cc.OF)
            
            r1cc.inst(r1ox.m_dot_ox)
            r1ox.inst(r1cc.P_cc)

            #add to arrays
            time_arr.append(r1ox.t)
            m_dot_arr.append(r1ox.m_dot_ox)
            thrust_arr.append(r1cc.instThrust)
            p_cc_arr.append(r1cc.P_cc) #NOTE: convert to atmospheric depending on data
            p_tank_arr.append(r1ox.P_tank)

            total_propellant = r1cc.total_propellant


        #print(r1ox.P_tank, r1cc.P_cc,r1ox.P_tank- r1cc.P_cc,  )

        #print(r1ox.m_ox, r1cc.m_fuel_t)

        #print(r1ox.t, r1cc.v_exit,r1cc.m_dot_cc_t,r1cc.R)

        #print("time: ", r1ox.t)

    print("\n", "### mass balance ###")
    print("total propellant calculated through sim: ", total_propellant+r1ox.m_ox+r1cc.m_fuel_t, "total starting/input propellant: ", inputs.m_ox+inputs.m_fuel_i, "difference (conservation of mass): ",total_propellant+r1ox.m_ox+r1cc.m_fuel_t -inputs.m_ox-inputs.m_fuel_i)

    ###WRITE CSV FOR FLIGHT SIM AND VALIDATION
    to_csv(time_arr,m_dot_arr, "m_dot_ox")
    to_csv(time_arr,thrust_arr, "thrust")
    to_csv(time_arr,p_cc_arr, "p_cc")
    to_csv(time_arr,p_tank_arr, "p_tank")


    ### LIQUID SETUP
    if len(inputs.analysis_mode) == 3:

        #input cc model
        cc_model_class = get_model(inputs.analysis_mode[0],'C')
        r1cc = cc_model_class(inputs.oxName, inputs.fuelName, inputs.CEA_fuel_str, inputs.m_fuel_i, 
                inputs.rho_fuel, inputs.a, inputs.n, inputs.L, inputs.A_port_i, 
                inputs.P_atm, inputs.A_throat, inputs.A_exit, inputs.timestep,)

        #input ox tank model
        OxTank_model_class = get_model(inputs.analysis_mode[1], 'T')
        r1ox = OxTank_model_class(inputs.oxName, inputs.timestep, inputs.m_ox, inputs.C_inj,
                inputs.V_tank, inputs.P_tank, inputs.P_atm, inputs.all_error, inputs.inj_model)

        #input fuel tank model
        FuelTank_model_class = get_model(inputs.analysis_mode[2], 'T')
        r1ox = OxTank_model_class(inputs.oxName, inputs.timestep, inputs.m_ox, inputs.C_inj,
                inputs.V_tank, inputs.P_tank, inputs.P_atm, inputs.all_error, inputs.inj_model)

    ### LIQUID THRUST CURVE

    ###PLOTS
    if inputs.thrust_curve_graphs == True:
        plt.subplot(1,4,1)
        plt.plot(time_arr,m_dot_arr)
        plt.xlabel('Time (s)')
        plt.ylabel('m_dot_ox (kg/s)')
        plt.title('Mass Flow Rate Over Time')
        plt.grid(True)

        plt.subplot(1,4,2)
        plt.plot(time_arr,thrust_arr)
        plt.xlabel('Time (s)')
        plt.ylabel('Thrust (N)')
        plt.title('Thrust Curve')
        plt.grid(True)

        plt.subplot(1,4,3)
        plt.plot(time_arr,p_cc_arr)
        plt.xlabel('Time (s)')
        plt.ylabel('Chamber Pressure (Pa)')
        plt.title('Chamber Pressure Over Time')
        plt.grid(True)

        plt.subplot(1,4,4)
        plt.plot(time_arr,p_tank_arr)
        plt.xlabel('Time (s)')
        plt.ylabel('Tank Pressure (Pa)')
        plt.title('Tank Pressure Over Time')
        plt.grid(True)


        plt.show()