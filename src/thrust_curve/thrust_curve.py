CC_MODEL_MAP = {
    1: 'src.models.hybrid_cc_w_fuel_grain',
    2: 'src.models.adiabatic_lre_cc'
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
import inspect

def get_variable_name(var):
    """Get the name of the variable as a string."""
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


    print(str,module_path)    
    if module_path is None:
        raise ValueError(f"Model {model_code} is not defined.")
    
    print("HERE!!!", module_path)
    module = importlib.import_module(module_path)
    return module.model  # or whichever class is appropriate for your tank model



def run_thrust_curve(inputs):

    ### SETUP FOR BOTH ENGINES:
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
                inputs.P_atm, inputs.A_throat, inputs.A_exit, inputs.timestep,)
        
    if inputs.analysis_mode[0] == 2:
        r1cc = cc_model_class(inputs.oxidizer_name, inputs.fuel_name, inputs.A_throat, inputs.A_exit, inputs.P_atm, inputs.TIMESTEP)


    ### Oxidizer Tank Setup:
    OxTank_model_class = get_model(inputs.analysis_mode[1], 'T')

    if inputs.analysis_mode[1] == 1:
        r1ox = OxTank_model_class(inputs.oxName, inputs.timestep, inputs.m_ox, inputs.C_inj,
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
        while (r1ox.t < 20):
            #print(r1cc.OF)
            r1cc.inst(r1ox.m_dot_ox)
            r1ox.inst(r1cc.P_cc)

            #add to arrays
            time_arr.append(r1ox.t)
            m_dot_arr.append(r1ox.m_dot_ox)
            thrust_arr.append(r1cc.instThrust)
            p_cc_arr.append(r1cc.P_cc) #NOTE: convert to atmospheric depending on data
            p_ox_tank_arr.append(r1ox.P_tank)
            total_propellant = r1cc.total_propellant


            #print(r1ox.P_tank, r1cc.P_cc,r1ox.P_tank- r1cc.P_cc,  )

            #print(r1ox.m_ox, r1cc.m_fuel_t)

            #print(r1ox.t, r1cc.v_exit,r1cc.m_dot_cc_t,r1cc.R)

            print("time: ", r1ox.t)

        print("\n", "### mass balance ###")
        print("total propellant calculated through sim: ", total_propellant+r1ox.m_ox+r1cc.m_fuel_t, "total starting/input propellant: ", inputs.m_ox+inputs.m_fuel_i, "difference (conservation of mass): ",total_propellant+r1ox.m_ox+r1cc.m_fuel_t -inputs.m_ox-inputs.m_fuel_i)


    ### for liquid setup fuel tank
    if len(inputs.analysis_mode) == 3:
        pressure_data = [p_cc_arr, p_ox_tank_arr, p_fuel_tank_arr]

        if inputs.analysis_mode[1] == 1:
            print("model invalid for fuel tank")
        
        if inputs.analysis_mode[1] == 2:
            print("model invalid for fuel tank")
            
        if inputs.analysis_mode[1] == 3:
            r1ox = OxTank_model_class(inputs.pressurant_name, inputs.m_pressurant, inputs.fuel_name, inputs.m_fuel, inputs.P_fueltank, inputs.ID_PROPTANK, inputs.TIMESTEP)
        
    ### LIQUID THRUST CURVE
    """
    time_arr = []
    m_dot_arr = []
    thrust_arr = []
    p_cc_arr = []
    p_ox_tank_arr = []
    p_fuel_tank_arr = []
    pressure_data = [p_cc_arr, p_tank_arr]
    P_cc = inputs.P_atm
    while(pressurantTank.P_prestank > pressurantTank.P_proptank):
    #while(t<5*TIMESTEP):

        #RECORD DATA
        time_arr.append(t)
        m_expelled += pressurantTank.m_dot*TIMESTEP
        m_pres_arr.append(pressurantTank.m_pres)
        P_tank_arr.append(pressurantTank.P_prestank)
        P_proptank_arr.append(pressurantTank.P_proptank)
        m_dot_arr.append(pressurantTank.m_dot)
        v_tank_arr.append(1/pressurantTank.rho_pres)
        T_tank_arr.append(pressurantTank.T_prestank)
        s_tank_arr.append(pressurantTank.s_prestank)

        #RUN THROUGH EACH CV AT EACH INSTANT
        pressurantTank.inst(oxidizerTank.P_proptank)
        oxidizerTank.inst(pressurantTank.m_dot_pres, pressurantTank.h_pres, t)

        #TODO: ADD FUEL TANK MODEL -->
        fuelTank.inst(P_atm)
        #TODO: ADD COMBUSTION CHAMBER MODEL -->`
        lre_cc.inst(1,2) 
        """


    ###WRITE CSV FOR FLIGHT SIM, VALIDATION AND OTHER EXTERNAL ANALYSIS
    to_csv(time_arr,m_dot_arr, "m_dot_ox")
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

        plt.show()

