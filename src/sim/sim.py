import csv
import importlib
import inspect
import numpy as np
import pandas as pd
from src.utils.model_registry import *

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


def get_model(kind: str, code: int):
    if kind == "T":
        module_path = OX_TANK_MODEL_MAP.get(code)
    elif kind == "F":
        module_path = FUEL_TANK_MODEL_MAP.get(code)
    elif kind == "I":
        module_path = INJECTOR_MODEL_MAP.get(code)
    elif kind == "C":
        module_path = CC_MODEL_MAP.get(code)
    elif kind == "N":
        module_path = NOZZLE_MODEL_MAP.get(code)
    else:
        raise ValueError(f"Unknown model kind {kind}")

    if module_path is None:
        raise ValueError(f"Model {code} is not defined in {kind} map")

    module_path, class_name = module_path.rsplit(".", 1)
    module = importlib.import_module(module_path)
    return getattr(module, class_name)



def check_nan(value):
    return 9e5 if np.isnan(value) else value




def sim(inputs):
    mode = inputs.mode
    timestep = inputs.TIMESTEP
    sim_time = inputs.sim_time

    # Build components if available
    ox_tank, fuel_tank, cc, nozzle = None, None, None, None

    if hasattr(inputs, "ox_tank_model"):
        TankClass = get_model("T", inputs.ox_tank_model)
        InjClass  = get_model("I", inputs.ox_inj_model)

        ox_inj    = InjClass(**inputs.ox_inj_kwargs)
        ox_tank   = TankClass(injector=ox_inj, timestep=timestep, **inputs.ox_tank_kwargs)

    if hasattr(inputs, "fuel_tank_model"):
        FuelTankClass = get_model("F", inputs.fuel_tank_model)
        fuel_tank = FuelTankClass(timestep=timestep, **inputs.fuel_tank_kwargs)
        # SPI injector is hardcoded inside fuel tank class

    if hasattr(inputs, "nozzle_model"):
            NozClass = get_model("N", inputs.nozzle_model)
            nozzle = NozClass(**inputs.nozzle_kwargs)

    if hasattr(inputs, "cc_model"):
        CCClass = get_model("C", inputs.cc_model)
        cc = CCClass(nozzle=nozzle, P_atm=inputs.P_atm, timestep=timestep, **inputs.cc_kwargs)

        # you don’t need to rebuild nozzle again here
        # nozzle already exists above if nozzle_model was defined

    t = 0.0
    records = []
    P_cc = inputs.P_atm

    if mode == "ox_tank":
        while t < sim_time:
            ox_out = ox_tank.inst(P_cc) or {}
            record = {"time": t, **ox_out}
            records.append(record)
            t += timestep

    elif mode == "fuel_tank":
        while t < sim_time:
            fuel_out = fuel_tank.inst(P_cc) or {}
            record = {"time": t, **fuel_out}
            records.append(record)
            t += timestep

    elif mode == "full_stack":
        # default dicts so they're always defined
        ox_out   = {"m_dot_ox": 0.0, "P_ox_tank": None}
        fuel_out = {"m_dot_fuel": 0.0, "P_fuel_tank": None}
        cc_out   = {"P_cc": P_cc, "thrust": 0.0, "m_dot_cc": 0.0}

        while t < sim_time:
            if ox_tank:
                ox_out = ox_tank.inst(P_cc) or ox_out
            if fuel_tank:
                fuel_out = fuel_tank.inst(P_cc) or fuel_out

            cc_out = cc.inst(
                ox_out["m_dot_ox"],
                fuel_out.get("m_dot_fuel", 0.0),
            ) or cc_out

            record = {"time": t, **ox_out, **fuel_out, **cc_out}
            records.append(record)

            P_cc = cc_out.get("P_cc", P_cc)  # feedback chamber pressure to tanks
            t += timestep

    df = pd.DataFrame(records)
    return df





















"""
def run_thrust_curve(inputs):

    ### SETUP FOR BOTH ENGINES:

    r1cc = None
    r1ox = None
    s1_fuel_tank = None

     t = 0.0
    time_arr, thrust_arr, P_cc_arr, P_ox_arr, P_fuel_arr, m_dot_arr = [], [], [], [], [], []
    
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
        while (r1ox.t < inputs.sim_time):
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
            s1_fuel_tank = fuel_tank_model_class(inputs.pressurant_name, inputs.m_pressurant, inputs.fuel_name, inputs.m_fuel, inputs.P_fueltank, inputs.ID_PROPTANK, inputs.V_tank_2, inputs.Cd_2, inputs.A_inj_2, inputs.T_amb, inputs.TIMESTEP)
        
        if inputs.analysis_mode[2] == 4:
            s1_fuel_tank = fuel_tank_model_class(inputs.T_tank, inputs.Cd_2, inputs.A_inj_2, inputs.fuel_name, inputs.fuel_tank_pressure_filepath, inputs.TIMESTEP)

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
            #print(r1cc.OF)
            p = r1cc.P_cc
            if r1cc.P_cc == 0 or r1cc.P_cc == -float('nan'):
                p = 9e4
            #print(p)
            r1ox.inst(p)
            s1_fuel_tank.inst(p)

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

                #p_fuel_up_min_dp = s1_fuel_tank.P_tank
                #p_fuel_down_min_dp = r1cc.P_cc

                #m_dot_fuel_min_dp = s1_fuel_tank.m_dot_fuel
                #rho_fuel_min_dp = s1_fuel_tank.rho_prop
                #kinematic_visc_fuel_min_dp = s1_fuel_tank.kinematic_visc_fuel
                #y_fuel_min_dp = s1_fuel_tank.y_fuel

                #t_fuel_min_dp = r1ox.t



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


        #print(f"\nThroat Properties at Peak Thrust for Heat Transfer\n------------\nRatio of specific heats: {y_peak} (-)\nSpec. Heat Const. Pres. {cp_peak} (J/(kg K))\nThroat Pressure {P_cc_peak} (Pa)\nCharacteristic Velocity {C_star_peak} (m/s)\nThroat Flame Temp {T_flame_peak} (K)\nViscosity {viscosity_peak} (Pa s)\nGas Constant {R_peak} (J/(kg K))")

        #print(f"\nMinimum Pressure Drop Fuel Inj Properties for Sizing\n------------\nTotal Fuel Mass Flow rate of all elements: {m_dot_fuel_min_dp} (kg/s)\nUpstream Pressure at inst: {p_fuel_up_min_dp} (Pa)\nDownstream Pressure at inst: {p_fuel_down_min_dp} (Pa)\nFuel Density at Orifice Outlet {rho_fuel_min_dp} (kg/m^3)\nFuel Kinematic Viscosity {kinematic_visc_fuel_min_dp} (Pa s)\nFuel Ratio of specific heats: {y_fuel_min_dp} (-)\nFuel Orifice Discharge Coeff: {inputs.Cd_2} (-)\nAt t = {t_fuel_min_dp} (s)")
        #print(f"\nMinimum Pressure Drop Ox Inj Properties for Sizing\n------------\nTotal Ox Mass Flow rate of all elements: {m_dot_ox_min_dp} (kg/s)\nUpstream Pressure at inst: {p_ox_up_min_dp} (Pa)\nDownstream Pressure at inst: {p_ox_down_min_dp} (Pa)\nOx Density at Orifice Outlet {rho_ox_min_dp} (kg/m^3)\nOx Kinematic Viscosity {kinematic_visc_ox_min_dp} (Pa s)\nOx Ratio of specific heats: {y_ox_min_dp} (-)\nOx Orifice Discharge Coeff: {inputs.Cd_1} (-)\nAt t = {t_ox_min_dp} (s)")

        plt.show()
"""

"""
from tqdm import tqdm  # add at the top

def sim(inputs):
    mode = inputs.mode
    timestep = inputs.TIMESTEP
    sim_time = inputs.sim_time

    # Build components if available
    ox_tank, fuel_tank, cc, nozzle = None, None, None, None
    ...
    t = 0.0
    records = []
    P_cc = inputs.P_atm

    steps = int(sim_time / timestep)

    if mode == "ox_tank":
        for _ in tqdm(range(steps), desc="Simulating ox tank"):
            ox_out = ox_tank.inst(P_cc) or {}
            record = {"time": t, **ox_out}
            records.append(record)
            t += timestep

    elif mode == "fuel_tank":
        for _ in tqdm(range(steps), desc="Simulating fuel tank"):
            fuel_out = fuel_tank.inst(P_cc) or {}
            record = {"time": t, **fuel_out}
            records.append(record)
            t += timestep

    elif mode == "full_stack":
        ox_out   = {"m_dot_ox": 0.0, "P_ox_tank": None}
        fuel_out = {"m_dot_fuel": 0.0, "P_fuel_tank": None}
        cc_out   = {"P_cc": P_cc, "thrust": 0.0, "m_dot_cc": 0.0}

        for _ in tqdm(range(steps), desc="Simulating full stack"):
            if ox_tank:
                ox_out = ox_tank.inst(P_cc) or ox_out
            if fuel_tank:
                fuel_out = fuel_tank.inst(P_cc) or fuel_out

            cc_out = cc.inst(
                ox_out["m_dot_ox"],
                fuel_out.get("m_dot_fuel", 0.0),
            ) or cc_out

            record = {"time": t, **ox_out, **fuel_out, **cc_out}
            records.append(record)

            P_cc = cc_out.get("P_cc", P_cc)
            t += timestep

    df = pd.DataFrame(records)
    return df


"""