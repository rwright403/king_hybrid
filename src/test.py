from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel
import numpy as np
import constants

add_new_fuel(constants.fuelName, constants.CEA_fuel_str)
C = CEA_Obj(oxName=constants.oxName, fuelName=constants.fuelName, pressure_units='Pa', isp_units='sec', cstar_units='m/s', temperature_units='K', sonic_velocity_units='m/s',enthalpy_units='kJ/kg',density_units='kg/m^3',specific_heat_units='kJ/kg-K')

expratio = constants.A_exit/constants.A_throat

P_cc = 4e6
CR = 0.01

while CR < 4:
    M_cc = C.get_Chamber_MachNumber(P_cc,6,CR) #TODO: ADD CR TO CONSTANTS!!!!
    print(CR, M_cc)
    CR+=0.1
    