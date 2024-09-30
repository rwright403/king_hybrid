CC_MODEL_MAP = {
    1: 'src.models.hybrid_cc_w_fuel_grain',
    2: 'src.models.adiabatic_lre_cc'
}

TANK_MODEL_MAP = {
    1: 'src.models.bens_ox_tank',
    2: 'src.models.adiabatic_ext_pressure_fed_cryo',
    3: 'src.models.adiabatic_pressurized_liquid_tank'
}


import importlib
   
import matplotlib.pyplot as plt
import numpy as np



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



class DataClass:
    def __init__(self, time_arr, m_dot_arr, thrust_arr, p_cc_arr, p_ox_tank_arr, p_fuel_tank_arr):
        self.time_arr_ = time_arr
        self.m_dot_arr_ = m_dot_arr
        self.thrust_arr_ = thrust_arr
        self.p_cc_arr_ = p_cc_arr
        self.p_ox_tank_arr_ = p_ox_tank_arr
        self.p_fuel_tank_arr_ = p_fuel_tank_arr


#TODO: start by calculating variance

def sensitivityAnalysis(inputs):
    
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
                inputs.P_atm, inputs.A_throat, inputs.A_exit, inputs.timestep,)
        
    if inputs.analysis_mode[0] == 2:
        r1cc = cc_model_class(inputs.oxidizer_name, inputs.fuel_name, inputs.A_throat, inputs.A_exit, inputs.P_atm, inputs.TIMESTEP)


    ### Oxidizer Tank Setup:
    OxTank_model_class = get_model(inputs.analysis_mode[1], 'T')

    if inputs.analysis_mode[1] == 1:
        r1ox = OxTank_model_class(inputs.oxName, inputs.timestep, inputs.m_ox, inputs.Cd_1, inputs.A_inj_1,
                inputs.V_tank, inputs.P_tank, inputs.P_atm, inputs.all_error, inputs.inj_model)
    
    if inputs.analysis_mode[1] == 2:
        pass
        #TODO: how to handle two functions!!!!
        #pressurantTank = simpleAdiabaticExtPressurantTank(pressurant_name, P_prestank, 0.5, P_oxtank, 0.01, OUTLET_DIAM,TIMESTEP)
        #oxidizerTank = simpleAdiabaticPressurizedTank(oxidizer_name, pressurant_name, ID_PROPTANK, P_oxtank, 5, V_PROPTANK, TIMESTEP)
    
    if inputs.analysis_mode[1] == 3:
        r1ox = OxTank_model_class(inputs.pressurant_name, inputs.m_pressurant, inputs.fuel_name, inputs.m_fuel, inputs.P_fueltank, inputs.ID_PROPTANK, inputs.V_tank_2, inputs.Cd_2, inputs.A_inj_2, inputs.TIMESTEP)
        
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


        iteration = DataClass(time_arr,m_dot_arr,thrust_arr,p_cc_arr,p_ox_tank_arr,None)

        return iteration


    ### for liquid setup fuel tank
    if len(inputs.analysis_mode) == 3:
        fuel_tank_model_class = get_model(inputs.analysis_mode[2],'T')

        pressure_data = [p_cc_arr, p_ox_tank_arr, p_fuel_tank_arr]

        if inputs.analysis_mode[2] == 1:
            print("model invalid for fuel tank (cannot use hybrid cc for liquid engine)")
        
        if inputs.analysis_mode[2] == 2:
            print("todo: implement")
            
        if inputs.analysis_mode[2] == 3:
            s1_fuel_tank = fuel_tank_model_class(inputs.pressurant_name, inputs.m_pressurant, inputs.fuel_name, inputs.m_fuel, inputs.P_fueltank, inputs.ID_PROPTANK, inputs.V_tank_2, inputs.Cd_2, inputs.A_inj_2, inputs.T_amb, inputs.TIMESTEP)
        
        ### LIQUID THRUST CURVE
        
        time_arr = []
        m_dot_arr = []
        thrust_arr = []
        p_cc_arr = []
        p_ox_tank_arr = []
        p_fuel_tank_arr = []
        pressure_data = [p_cc_arr, p_ox_tank_arr, p_fuel_tank_arr]
        P_cc = inputs.P_atm

        #SAVE VALUES TO PRINT OUTPUTS FOR HEAT TRANSFER HAND CALC:
        """
        y_peak = 0
        cp_peak = 0
        P_cc_peak = 0
        C_star_peak = 0
        T_flame_peak = 0
        viscosity_peak = 0
        """

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
    
            m_fuel_burned += s1_fuel_tank.m_dot_fuel*r1ox.timestep
            m_ox_burned += r1ox.m_dot_ox*r1ox.timestep

            if (r1cc.instThrust > r1cc.prev_thrust):

                
                arr1 = r1cc.C.get_Throat_MolWt_gamma(r1cc.P_cc, r1cc.OF, r1cc.expratio, frozen=0)
                arr2 = r1cc.C.get_Temperatures(r1cc.P_cc, r1cc.OF, r1cc.expratio, frozen=0, frozenAtThroat=0)
                arr3 = r1cc.C.get_Throat_Transport(r1cc.P_cc, r1cc.OF, r1cc.expratio, frozen=0)
                
                """
                y_peak = arr1[1]
                cp_peak = arr3[0]
                P_cc_peak = r1cc.P_cc*(2/(y_peak+1))**(y_peak/(y_peak-1))
                C_star_peak = r1cc.C.get_Cstar(r1cc.P_cc, r1cc.OF)
                T_flame_peak = arr2[1]
                viscosity_peak = arr3[1]
                R_peak = r1cc.R
                """
            #add to arrays
            time_arr.append(r1ox.t)
            m_dot_arr.append(r1ox.m_dot_ox)
            thrust_arr.append(r1cc.instThrust)
            p_cc_arr.append(r1cc.P_cc) #NOTE: convert to atmospheric depending on data
            p_ox_tank_arr.append(r1ox.P_tank)
            p_fuel_tank_arr.append(s1_fuel_tank.P_tank)

        iteration = DataClass(time_arr,m_dot_arr,thrust_arr,p_cc_arr,p_ox_tank_arr,p_fuel_tank_arr)
        #print(f"\n{len(time_arr)}\n")
        return iteration
    
    

def produce_graphs(big_data,i_arr,It_arr):
    j=1
    while(j<6):

        if(j==1):
            k=0
            for s in big_data:
                plt.subplot(1,5,j)
                plt.plot(s.time_arr_,s.m_dot_arr_, label= f'{i_arr[k]:.5}')
                k+=1
            plt.legend()
            plt.xlabel('Time (s)')
            plt.ylabel('m_dot_ox (kg/s)')
            plt.title('Mass Flow Rate Over Time')
            plt.grid(True)

        if(j==2):
            k=0
            for s in big_data:
                plt.subplot(1,5,j)
                plt.plot(s.time_arr_,s.thrust_arr_, label= f"Total Impulse: {It_arr[k]:.6}")
                k+=1
            plt.legend()
            plt.xlabel('Time (s)')
            plt.ylabel('Thrust (N)')
            plt.title('Thrust Curve')
            plt.grid(True)

        if(j==3):
            k=0
            for s in big_data:
                plt.subplot(1,5,j)
                plt.plot(s.time_arr_,s.p_cc_arr_)
                k+=1
            plt.xlabel('Time (s)')
            plt.ylabel('Chamber Pressure (Pa)')
            plt.title('Chamber Pressure Over Time')
            plt.grid(True)

        if(j==4):
            k=0
            for s in big_data:
                plt.subplot(1,5,j)
                plt.plot(s.time_arr_,s.p_ox_tank_arr_)
                k+=1
            plt.xlabel('Time (s)')
            plt.ylabel('Oxidizer Tank Pressure (Pa)')
            plt.title('Oxidizer Tank Pressure Over Time')
            plt.grid(True)

        ###TODO: add logic so this only appears if simulating a lre
        if(j==5):
            k=0
            for s in big_data:
                plt.subplot(1,5,j)
                plt.plot(s.time_arr_,s.p_fuel_tank_arr_)
                k+=1
            plt.xlabel('Time (s)')
            plt.ylabel('Fuel Tank Pressure (Pa)')
            plt.title('Fuel Tank Pressure Over Time')
            plt.grid(True)
            

        j += 1
    plt.show()
    


def update_i(inputs):
    return (getattr(inputs,inputs.test_var_name) + (inputs.max_bound-inputs.min_bound)/(inputs.num_iterations -1) )

def run_sensitivity_analysis(inputs):

    big_data = []
    i_arr = []
    It_arr = []

    if hasattr(inputs, inputs.test_var_name):

        #update in the module dynamically using settart()
        setattr(inputs, inputs.test_var_name, inputs.min_bound)

        while(getattr(inputs,inputs.test_var_name)<=inputs.max_bound):
            
            results = sensitivityAnalysis(inputs)
            big_data.append( results )

            total_impulse = np.trapz(results.thrust_arr_, results.time_arr_)
            It_arr.append(total_impulse)

            i_arr.append(getattr(inputs,inputs.test_var_name) + (inputs.max_bound-inputs.min_bound)/(inputs.num_iterations -1) )
            
            #update variable
            i = (getattr(inputs,inputs.test_var_name) + (inputs.max_bound-inputs.min_bound)/(inputs.num_iterations -1) )
            setattr(inputs,inputs.test_var_name, i)

        #print(f"\n{len(big_data)}\n")
            
        produce_graphs(big_data,i_arr,It_arr)

        #NOTE: what is i_arr used for again?    
    
















    #malding if statements, you hate to see it for each inputs, then call function with sensitivity analysis inside
"""
    #TODO: PROBABLY NICE TO PRINT OUT CONSTANTS
        
    print("\nSensitivity Analysis Summary:\n")
    print(f"test_var = \"{inputs.test_var}\"")
    print(f"min_bound = {inputs.min_bound}")
    print(f"max_bound = {inputs.max_bound}")
    print(f"num_iterations = {inputs.num_iterations}\n")
    print(f"ENGINE DATA\n")
    print(f"oxName = '{inputs.oxName}'")
    print(f"fuelName = '{inputs.fuelName}'")
    print(f"rho_fuel = {inputs.rho_fuel}")
    print(f"m_fuel_i = {inputs.m_fuel_i}")
    print(f"a = {inputs.a}")
    print(f"n = {inputs.n}")
    print(f"L = {inputs.L}")
    print(f"A_port_i = {inputs.A_port_i}")
    print(f"A_throat = {inputs.A_throat}")
    print(f"A_exit = {inputs.A_exit}")
    print(f"r_tank = {inputs.r_tank}")
    print(f"height_tank = {inputs.height_tank}")
    print(f"V_tank = {inputs.V_tank}")
    print(f"P_tank = {inputs.P_tank}")
    print(f"m_ox = {inputs.m_ox}")
    print(f"C_inj = {inputs.C_inj}")
    print(f"P_atm = {inputs.P_atm}")
    print(f"timestep = {inputs.timestep}")
    print(f"sim_time = {inputs.sim_time}")
    #i love good code :)

    plt.show()
"""