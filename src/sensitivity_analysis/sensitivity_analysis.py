###TODO: ADD PCC AS INPUT!!!!!!
from src import constants

from src.thrust_curve.bens_ox_tank import OxTank
from src.thrust_curve.combustion_chamber import cc
   
import matplotlib.pyplot as plt
import numpy as np


class DataClass:
    def __init__(self, time_arr, m_dot_arr, thrust_arr, p_cc_arr, p_tank_arr):
        self.time_arr = time_arr
        self.m_dot_arr = m_dot_arr
        self.thrust_arr = thrust_arr
        self.p_cc_arr = p_cc_arr
        self.p_tank_arr = p_tank_arr




#TODO: start by calculating variance

def sensitivityAnalysis(fill_level,C_inj,V_tank,P_tank,m_fuel_i, a, n, L, A_port_i, A_throat, A_exit, P_cc):

    i = constants.min_bound

    time_arr = []
    m_dot_arr = []
    thrust_arr = []
    p_cc_arr = []
    p_tank_arr = []

    big_data = []

    while(i<constants.max_bound):
        #classes give instantaneous data for giga struct
        print("in sens while loop")

        r1ox = OxTank(constants.oxName, constants.timestep, fill_level, C_inj, V_tank, P_tank, constants.P_atm, constants.all_error)

        r1cc = cc(constants.oxName, constants.fuelName, constants.CEA_fuel_str, m_fuel_i, 
            constants.rho_fuel, a, n, L, A_port_i,constants.P_atm, A_throat, A_exit, constants.timestep)
        
        while r1ox.t < 7:
            r1ox.inst(P_cc)
            r1cc.inst(r1ox.m_dot_ox, P_cc)

            time_arr.append(r1ox.t)
            m_dot_arr.append(r1ox.m_dot_ox)
            thrust_arr.append(r1cc.instThrust)
            p_cc_arr.append(r1cc.P_cc)
            p_tank_arr.append(r1ox.P_tank)
        
        #probably add arrays to bitchass array of classes
        iteration = DataClass(time_arr,m_dot_arr,thrust_arr,p_cc_arr,p_tank_arr)
        big_data.append(iteration)
        i = i + (constants.max_bound-constants.min_bound)/constants.num_iterations

    print("done")


#malding if statements, you hate to see it for each input, then call function with sensitivity analysis inside

if constants.test_var=="fill_level":
    sensitivityAnalysis(constants.min_bound,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit)

if constants.test_var=="C_inj":
    sensitivityAnalysis(constants.fill_level,constants.min_bound,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit)

if constants.test_var=="V_tank":
    sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.min_bound,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit)

if constants.test_var=="P_tank":
    sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.min_bound,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit)

if constants.test_var=="m_fuel_i":
    sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.min_bound,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit)

if constants.test_var=="a":
    sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.min_bound,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit)

if constants.test_var=="n":
    sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.min_bound,constants.L,constants.A_port_i, constants.A_throat, constants.A_exit)

if constants.test_var=="L":
    sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.min_bound,constants.A_port_i, constants.A_throat, constants.A_exit)

if constants.test_var=="A_port_i":
    sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.min_bound, constants.A_throat, constants.A_exit)

if constants.test_var=="A_throat":
    sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, constants.min_bound, constants.A_exit)

if constants.test_var=="A_exit":
    sensitivityAnalysis(constants.fill_level,constants.C_inj,constants.V_tank,constants.P_tank,constants.m_fuel_i,constants.a,constants.n,constants.L,constants.A_port_i, constants.A_throat, constants.min_bound)