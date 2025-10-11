import importlib

# Registry maps code --> import path
OX_TANK_MODEL_MAP = {
    1: "src.models.ox_tank.equilibrium_tank.equilibrium_tank_model",
    2: "src.models.ox_tank.non_equilibrium_tank.non_equilibrium_tank_model",
}


FUEL_TANK_MODEL_MAP = {
    1: "src.models.fuel_tank.pressurized_liq_fuel_tank.pressurized_liq_fuel_tank"
}

INJECTOR_MODEL_MAP = {
    1: "src.models.inj.spi.spi_model",
    2: "src.models.inj.hem.hem_model",
    3: "src.models.inj.dyer.dyer_model",
    4: "src.models.inj.modified_omega.modified_omega_model",
    5: "src.models.inj.dyer_hem_choke_pt.dyer_model"
}

CC_MODEL_MAP = {
    1: "src.models.cc.hybrid_cc_w_fuel_grain.hybrid_cc_w_fuel_grain_model",
    2: "src.models.cc.adiabatic_lre_cc.adiabatic_lre_cc_model",
    3: "src.models.cc.pipe_downstream_p_cc.pipe_downstream_p_cc"
}

NOZZLE_MODEL_MAP = {
    1: "src.models.nozzle.basic_nozzle.basic_nozzle_model",
}

DRAG_MODEL_MAP = {
    1: "",
    2: "src.models.drag.barrowman.barrowman",
}

def _load_class(path: str):
    #Dynamically load a class from a 'module.ClassName' string.
    module_path, class_name = path.rsplit(".", 1)
    module = importlib.import_module(module_path)
    return getattr(module, class_name)


def get_model(kind: str, code: int):
    #Return the class constructor for a model type (T/I/C/N/F).
    if kind == "T":
        return _load_class(OX_TANK_MODEL_MAP[code])
    if kind == "F":
        return _load_class(FUEL_TANK_MODEL_MAP[code])
    if kind == "I":
        return _load_class(INJECTOR_MODEL_MAP[code])
    if kind == "C":
        return _load_class(CC_MODEL_MAP[code])
    if kind == "N":
        return _load_class(NOZZLE_MODEL_MAP[code])
    raise ValueError(f"Unknown model kind {kind}")