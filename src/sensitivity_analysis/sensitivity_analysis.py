CC_MODEL_MAP = {
    1: 'src.models.hybrid_cc_w_fuel_grain',
    2: 'adiabatic_lre_cc'
}

TANK_MODEL_MAP = {
    1: 'src.models.bens_ox_tank',
    2: 'src.models.adiabatic_ext_pressure_fed_cryo',
    3: 'src.tank_models.adiabatic_pressurized'
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
    def __init__(self, time_arr, m_dot_arr, thrust_arr, p_cc_arr, p_tank_arr):
        self.time_arr_ = time_arr
        self.m_dot_arr_ = m_dot_arr
        self.thrust_arr_ = thrust_arr
        self.p_cc_arr_ = p_cc_arr
        self.p_tank_arr_ = p_tank_arr


#TODO: start by calculating variance

def sensitivityAnalysis(inputs):
    ### SETUP FOR BOTH ENGINES:

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
        

    ### for liquid setup fuel tank
    if len(inputs.analysis_mode) == 3:

        if inputs.analysis_mode[1] == 1:
            print("model invalid for fuel tank")
        
        if inputs.analysis_mode[1] == 2:
            print("model invalid for fuel tank")
            
        if inputs.analysis_mode[1] == 3:
            r1ox = OxTank_model_class(inputs.pressurant_name, inputs.m_pressurant, inputs.fuel_name, inputs.m_fuel, inputs.P_fueltank, inputs.ID_PROPTANK, inputs.TIMESTEP)
        

    time_arr = []
    m_dot_arr = []
    thrust_arr = []
    p_cc_arr = []
    p_tank_arr = []

    #P_cc = inputs.P_atm

    r1ox = OxTank(inputs.oxName, inputs.timestep, m_ox, C_inj, V_tank, P_tank, inputs.P_atm, inputs.all_error, inputs.inj_model)

    r1cc = cc(inputs.oxName, inputs.fuelName, inputs.CEA_fuel_str, inputs.m_fuel_i, 
        inputs.rho_fuel, a, n, L, A_port_i,inputs.P_atm, A_throat, A_exit, inputs.timestep)
        
    while r1ox.t < inputs.sim_time:
        r1ox.inst(r1cc.P_cc)
        r1cc.inst(r1ox.m_dot_ox)

        time_arr.append(r1ox.t)
        m_dot_arr.append(r1ox.m_dot_ox)
        thrust_arr.append(r1cc.instThrust) #this is returning 0 every time?? its because self.inst thrust = 0!
        p_cc_arr.append(r1cc.P_cc)
        p_tank_arr.append(r1ox.P_tank)

        iteration = DataClass(time_arr,m_dot_arr,thrust_arr,p_cc_arr,p_tank_arr)

    return iteration

def produce_graphs(big_data,i_arr):
    j=1
    while(j<5):

        if(j==1):
            k=0
            for s in big_data:
                plt.subplot(1,4,j)
                plt.plot(s.time_arr_,s.m_dot_arr_, label= f'{i_arr[k]}')
                k+=1
            plt.legend()
            plt.xlabel('Time (s)')
            plt.ylabel('m_dot_ox (kg/s)')
            plt.title('Mass Flow Rate Over Time')
            plt.grid(True)

        if(j==2):
            k=0
            for s in big_data:
                plt.subplot(1,4,j)
                plt.plot(s.time_arr_,s.thrust_arr_, label= f'{i_arr[k]}')
                k+=1
            plt.legend()
            plt.xlabel('Time (s)')
            plt.ylabel('Thrust (N)')
            plt.title('Thrust Curve')
            plt.grid(True)

        if(j==3):
            k=0
            for s in big_data:
                plt.subplot(1,4,j)
                plt.plot(s.time_arr_,s.p_cc_arr_, label= f'{i_arr[k]}')
                k+=1
            plt.legend()
            plt.xlabel('Time (s)')
            plt.ylabel('Chamber Pressure (Pa)')
            plt.title('Chamber Pressure Over Time')
            plt.grid(True)

        if(j==4):
            k=0
            for s in big_data:
                plt.subplot(1,4,j)
                plt.plot(s.time_arr_,s.p_tank_arr_, label= f'{i_arr[k]}')
                k+=1
            plt.legend()
            plt.xlabel('Time (s)')
            plt.ylabel('Tank Pressure (Pa)')
            plt.title('Tank Pressure Over Time')
            plt.grid(True)
            

        j += 1


def update_i(i):
    i_arr.append(i)
    return i + (inputs.max_bound-inputs.min_bound)/(inputs.num_iterations -1)


i = inputs.min_bound
big_data = []
i_arr = []
#malding if statements, you hate to see it for each inputs, then call function with sensitivity analysis inside

if inputs.test_var=="m_ox":
    while(i<=inputs.max_bound):
        big_data.append( sensitivityAnalysis(i,inputs.C_inj,inputs.V_tank,inputs.P_tank,inputs.m_fuel_i,inputs.a,inputs.n,inputs.L,inputs.A_port_i, inputs.A_throat, inputs.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if inputs.test_var=="C_inj":
    while(i<=inputs.max_bound):
        big_data.append( sensitivityAnalysis(inputs.m_ox,i,inputs.V_tank,inputs.P_tank,inputs.m_fuel_i,inputs.a,inputs.n,inputs.L,inputs.A_port_i, inputs.A_throat, inputs.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if inputs.test_var=="V_tank":
    while(i<=inputs.max_bound):
        big_data.append( sensitivityAnalysis(inputs.m_ox,inputs.C_inj,i,inputs.P_tank,inputs.m_fuel_i,inputs.a,inputs.n,inputs.L,inputs.A_port_i, inputs.A_throat, inputs.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if inputs.test_var=="P_tank":
    while(i<=inputs.max_bound):
        big_data.append( sensitivityAnalysis(inputs.m_ox,inputs.C_inj,inputs.V_tank,i,inputs.m_fuel_i,inputs.a,inputs.n,inputs.L,inputs.A_port_i, inputs.A_throat, inputs.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if inputs.test_var=="m_fuel_i":
    while(i<=inputs.max_bound):
        big_data.append( sensitivityAnalysis(inputs.m_ox,inputs.C_inj,inputs.V_tank,inputs.P_tank,i,inputs.a,inputs.n,inputs.L,inputs.A_port_i, inputs.A_throat, inputs.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if inputs.test_var=="a":
    while(i<=inputs.max_bound):
        big_data.append( sensitivityAnalysis(inputs.m_ox,inputs.C_inj,inputs.V_tank,inputs.P_tank,inputs.m_fuel_i,i,inputs.n,inputs.L,inputs.A_port_i, inputs.A_throat, inputs.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if inputs.test_var=="n":
    while(i<=inputs.max_bound):
        big_data.append( sensitivityAnalysis(inputs.m_ox,inputs.C_inj,inputs.V_tank,inputs.P_tank,inputs.m_fuel_i,inputs.a,i,inputs.L,inputs.A_port_i, inputs.A_throat, inputs.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if inputs.test_var=="L":
    while(i<=inputs.max_bound):
        big_data.append( sensitivityAnalysis(inputs.m_ox,inputs.C_inj,inputs.V_tank,inputs.P_tank,inputs.m_fuel_i,inputs.a,inputs.n,i,inputs.A_port_i, inputs.A_throat, inputs.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if inputs.test_var=="A_port_i":
    while(i<=inputs.max_bound):
        big_data.append( sensitivityAnalysis(inputs.m_ox,inputs.C_inj,inputs.V_tank,inputs.P_tank,inputs.m_fuel_i,inputs.a,inputs.n,inputs.L,i, inputs.A_throat, inputs.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if inputs.test_var=="A_throat":
    while(i<=inputs.max_bound):
        big_data.append( sensitivityAnalysis(inputs.m_ox,inputs.C_inj,inputs.V_tank,inputs.P_tank,inputs.m_fuel_i,inputs.a,inputs.n,inputs.L,inputs.A_port_i, i, inputs.A_exit) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)

if inputs.test_var=="A_exit":
    while(i<=inputs.max_bound):
        big_data.append( sensitivityAnalysis(inputs.m_ox,inputs.C_inj,inputs.V_tank,inputs.P_tank,inputs.m_fuel_i,inputs.a,inputs.n,inputs.L,inputs.A_port_i, inputs.A_throat, i) )
        i = update_i(i)
    produce_graphs(big_data,i_arr)


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