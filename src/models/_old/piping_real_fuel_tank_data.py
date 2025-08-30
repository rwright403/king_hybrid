import CoolProp.CoolProp as CP #I love coolprop! ~ units: http://www.coolprop.org/v4/apidoc/CoolProp.html
import matplotlib.pyplot as plt
import numpy as np
import csv

def read_csv(file_path):   
    x = []
    y = []
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:

            if not row:
                continue

            try:
                x.append(float(row[0]))
                y.append(float(row[1]))
            except ValueError:
                continue

    return np.array(x), np.array(y)

class model():
    def __init__(self, T_tank, Cd_spi, A_inj, fuel_name, fuel_tank_pressure_filepath,TIMESTEP): #want tank temperature for coolprop

        self.T_tank = T_tank
        self.A_inj = A_inj
        self.Cd_spi = Cd_spi
        self.fuel_name = fuel_name
        self.P_tank = 0
        self.t = 0
        self.m_dot_fuel = 0
        self.TIMTESTEP = TIMESTEP

        #copy data from csv in inputs.py file and store here
        self.Time_arr, self.Tank_Pres_arr = read_csv(fuel_tank_pressure_filepath)


    def inst(self,P_cc): #want downstream pressure and time?

        #interpolate data from vectors to get real time pressure for spi model
        self.P_tank = np.interp(self.t, self.Time_arr, self.Tank_Pres_arr)

        #call coolprop to solve density at inlet for spi model
        #NOTE: this works assuming piston fed systems the nos is a subcooled liquid which i think is a good assumption
        rho_tank_exit = CP.PropsSI('D', 'T', self.T_tank, 'P', self.P_tank, self.fuel_name)

        #return mass flow rate with spi model:
        if self.P_tank <= P_cc:
            self.m_dot_fuel = 0
        else:
            self.m_dot_fuel = self.Cd_spi * self.A_inj * np.sqrt( 2 * rho_tank_exit * (self.P_tank - P_cc)  )
        self.t += self.TIMTESTEP
        #print(self.m_dot_fuel)
        return self.m_dot_fuel