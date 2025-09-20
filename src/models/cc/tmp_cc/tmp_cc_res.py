import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel
from src.models.nozzle._base import BaseNozzle


def area_mach_func(M, gamma, A_ratio):
    term1 = 1 / M
    term2 = (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
    term3 = (1 + (gamma-1)/2 * M**2)**((gamma+1)/(2*(gamma-1)))
    return (term1 * term2 * term3) - A_ratio


class basic_nozzle_model(BaseNozzle):
    def __init__(self, d_throat: float, expratio: float, P_atm: float, C):
        super().__init__(d_throat, expratio, P_atm, C)

    def sol_thrust(self, P_cc, T_cc, gamma, R):
        Critical_PR = (2/(gamma+1)) ** (-gamma/(gamma+1))
        Actual_PR = self.P_atm / P_cc

        # Choked case
        if Actual_PR <= Critical_PR:
            M_exit = brentq(area_mach_func, 1.0001, 6,
                            args=(gamma, self.A_exit/self.A_throat))
            m_dot_exit = (self.A_throat * P_cc * np.sqrt(gamma/(R*T_cc)) *
                          ((gamma+1)/2) ** (-(gamma+1)/(2*(gamma-1))))
        else:
            # Subsonic case
            M_exit = brentq(area_mach_func, 1e-6, 0.99999,
                            args=(gamma, self.A_exit/self.A_throat))
            m_dot_exit = (self.A_exit * P_cc * np.sqrt(gamma/(R*T_cc)) *
                          (1 + (gamma-1)/2 * M_exit**2)**(-(gamma+1)/(2*(gamma-1))))

        T_exit = T_cc / (1 + (gamma-1)/2 * M_exit**2)
        a_exit = np.sqrt(R * gamma * T_exit)
        v_exit = a_exit * M_exit
        P_exit = P_cc / ((1 + (gamma-1)/2 * M_exit**2)**(gamma/(gamma-1)))

        self.instThrust = m_dot_exit*v_exit + self.A_exit*(P_exit - self.P_atm)
        return self.instThrust, m_dot_exit


def residual_curve(m_dot_ox, r, L, rho_f, a, n, nozzle, C, expratio):
    A_port = np.pi * r**2
    G_ox = m_dot_ox / A_port

    print("G_ox: ", G_ox)

    r_dot = a * G_ox**n
    m_dot_fuel = rho_f * (2*np.pi*r*L) * r_dot

    print("r_dot, m_dot_fuel: ", r_dot, m_dot_fuel)


    OF = m_dot_ox / m_dot_fuel

    print(OF, m_dot_ox, m_dot_fuel)

    pressures = np.linspace(5e5, 1e8, 200)
    residuals = []

    for Pc in pressures:

        print("cea inputs: ", Pc, OF, expratio)

        T_c = C.get_Tcomb(Pc, OF)
        MW, gamma = C.get_Chamber_MolWt_gamma(Pc, OF, expratio)
        R_spec = 8314.0 / MW

        print("noz inputs: ", Pc, T_c, gamma, R_spec)

        _, m_dot_noz = nozzle.sol_thrust(Pc, T_c, gamma, R_spec)
        residuals.append((m_dot_ox + m_dot_fuel) - m_dot_noz)

        print("pressure sweep: ",m_dot_noz, ((m_dot_ox + m_dot_fuel) - m_dot_noz),Pc  )

    return pressures, np.array(residuals)


# -----------------------
# Example usage
# -----------------------

# CEA object (use your oxidizer/fuel names here)
fuel_name = 'paraffin'
#C32H66 from RPA Paraffin Wax Composition
fuel_properties = f"""
fuel paraffin  C 32   H 66    wt%=100.00
h,KJ/Kgmol=-1860600     t(k)=298.15   rho,kg/m3={900}
"""

add_new_fuel(fuel_name, fuel_properties)
my_CEA = CEA_Obj(oxName="N2O", fuelName="paraffin")

# Basic nozzle
d_throat = 0.06   # 20 mm
expratio = 5.0
P_atm = 5e5
my_nozzle = basic_nozzle_model(d_throat, expratio, P_atm, my_CEA)

# Residual curve
pressures, residuals = residual_curve(
    m_dot_ox=1,    # kg/s
    r=0.03,          # m
    L=0.3,           # m
    rho_f=900,       # kg/m^3
    a=0.155e-3 , n=0.5,
    nozzle=my_nozzle,
    C=my_CEA,
    expratio=expratio
)

plt.semilogx(pressures, residuals)
plt.axhline(0, color='k', linestyle='--')
plt.xlabel("Chamber Pressure Pc [Pa]")
plt.ylabel("Mass balance residual [kg/s]")
plt.title("Inflow vs Nozzle Outflow Residual Curve")
plt.show()
