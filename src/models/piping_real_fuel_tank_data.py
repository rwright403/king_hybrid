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
    def __init__(self, T_tank, Cd_spi, A_inj, fuel_name, fuel_tank_pressure_filepath): #want tank temperature for coolprop

        self.T_tank = T_tank
        self.A_inj = A_inj
        self.Cd_spi = Cd_spi
        self.fuel_name = fuel_name

        #copy data from csv in inputs.py file and store here
        self.Time_arr, self.Tank_Pres_arr = read_csv(fuel_tank_pressure_filepath)


    def inst(self,P_cc, t): #want downstream pressure and time?

        #interpolate data from vectors to get real time pressure for spi model
        P_interpolated = np.interp(t, self.Time_arr, self.Tank_Pres_arr)

        #call coolprop to solve density at inlet for spi model
        #NOTE: this works assuming piston fed systems the nos is a subcooled liquid which i think is a good assumption
        rho_tank_exit = CP.PropsSI('D', 'T', self.T_tank, 'P', P_interpolated, self.fuel_name)

        #return mass flow rate with spi model:
        m_dot_spi = self.Cd_spi * self.A_inj * np.sqrt( 2 * rho_tank_exit * (P_interpolated - P_cc)  )
        return m_dot_spi