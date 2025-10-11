import csv
import importlib
import inspect
import numpy as np
import pandas as pd
from src.utils.model_registry import *
import traceback

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




def prop_sim(kwargs: dict):
    # ------------------------
    # Globals
    # ------------------------
    g = kwargs["globals_kwargs"]
    mode = g["mode"]
    timestep = g["timestep"]
    sim_time = g["sim_time"]
    P_atm = g["P_atm"]

    m = kwargs["models_kwargs"]

    # ------------------------
    # Build components
    # ------------------------
    ox_tank = None
    fuel_tank = None
    nozzle = None
    cc = None

    if kwargs.get("ox_tank_kwargs") is not None:
        TankClass = get_model("T", m["ox_tank_model"])
        InjClass  = get_model("I", m["ox_inj_model"])
        ox_inj    = InjClass(**kwargs["ox_inj_kwargs"])
        ox_tank   = TankClass(injector=ox_inj, timestep=timestep, **kwargs["ox_tank_kwargs"])

    if kwargs.get("fuel_tank_kwargs") is not None:
        FuelTankClass = get_model("F", m["fuel_tank_model"])
        fuel_tank = FuelTankClass(timestep=timestep, **kwargs["fuel_tank_kwargs"])

    if kwargs.get("nozzle_kwargs") is not None:
        NozClass = get_model("N", m["nozzle_model"])
        nozzle = NozClass(**kwargs["nozzle_kwargs"])


    cc = None

    if m["cc_model"] == 3:  # validation CC (needs P_cc file)
        if kwargs.get("cc_kwargs") is not None:
            CCClass = get_model("C", 3)
            cc = CCClass(**kwargs["cc_kwargs"])
        else:
            # no P_cc file --> skip CC instantiation
            cc = None

    else:  # all other CC models
        CCClass = get_model("C", m["cc_model"])
        base_kwargs = dict(nozzle=nozzle, P_atm=P_atm, timestep=timestep)
        cc = CCClass(**base_kwargs, **kwargs["cc_kwargs"])




    # ------------------------
    # Simulation loop
    # ------------------------
    t = 0.0
    records = []
    P_cc = P_atm

    try:

        if mode == "ox_tank":
            while t < sim_time:
                ox_out = ox_tank.inst(P_cc) or {}
                records.append({"time": t, **ox_out})

                if cc is not None:  # only update if a chamber exists
                    cc_out = cc.inst() or cc_out
                    P_cc = cc_out.get("P_cc", P_cc)
                    #print("P_cc: ", P_cc)
                # else: keep using the last P_cc (initialized to P_atm)

                t += timestep


        elif mode == "fuel_tank":
            while t < sim_time:
                fuel_out = fuel_tank.inst(P_cc) or {}
                records.append({"time": t, **fuel_out})
                t += timestep

        elif mode == "full_stack":

            ox_out   = {"m_dot_ox": 0.0, "P_ox_tank": None}
            fuel_out = {"m_dot_fuel": 0.0, "P_fuel_tank": None}
            cc_out   = {"P_cc": P_cc, "thrust": 0.0, "m_dot_fuel": 0.0}

            while t < sim_time:
                if ox_tank:
                    ox_out = ox_tank.inst(P_cc) or ox_out
                if fuel_tank:
                    fuel_out = fuel_tank.inst(P_cc) or fuel_out

                if cc:
                    cc_out = cc.inst(ox_out.get("m_dot_ox", 0.0),fuel_out.get("m_dot_fuel", 0.0)) #or cc_out

                records.append({"time": t, **ox_out, **fuel_out, **cc_out})
                P_cc = cc_out.get("P_cc", P_cc)

                t += timestep
        
        
    except Exception as e:
        print("Simulation error:", e)
        traceback.print_exc()

    # ------------------------
    # Cleanup / Resize arrays
    # ------------------------
    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)

    # Find minimum array length across all array-like cells
    min_len = None
    for col in df.columns:
        for val in df[col]:
            if hasattr(val, "__len__") and not isinstance(val, (str, bytes)):
                try:
                    l = len(val)
                    if min_len is None or l < min_len:
                        min_len = l
                except Exception:
                    pass
    if min_len is None:
        return df  # no array-like values, nothing to fix

    # Truncate array-like cells so everything is the same size
    for col in df.columns:
        def _fix(x):
            if hasattr(x, "__len__") and not isinstance(x, (str, bytes)):
                try:
                    return np.asarray(x)[:min_len]
                except Exception:
                    return x
            return x
        df[col] = df[col].apply(_fix)

    return df
