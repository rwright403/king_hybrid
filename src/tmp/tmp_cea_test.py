import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel

fuel_name       = "paraffin"
fuel_properties = """fuel paraffin  C 20   H 42    wt%=100.00
h,KJ/mol=-13313.49  t(k)=298.15   rho,kg/m3=900"""

fuel_properties = """
#fuel PARAFFIN_CH2  (CH2)x(s)  Ethylene polymer. Est. from paraffin HC series, TRC p-2500, 4/88.
C 1 H 2   wt%=100.00    h,cal=-25600.0    t(k)=298.15    rho=0.900"""

#add_new_fuel(fuel_name, fuel_properties)

# Initialize CEA object with your propellants
C = CEA_Obj(
    oxName="N2O",
    fuelName="paraffin",
    pressure_units="Pa",
    isp_units="sec",
    cstar_units="m/s",
    temperature_units="K",
    sonic_velocity_units="m/s",
    enthalpy_units="kJ/kg",
    density_units="kg/m^3",
    specific_heat_units="kJ/kg-K",
)

# Sweep some random chamber pressures and mixture ratios
pressures = [1e5, 5e5, 1e6, 2e6, 5e6]   # Pa
mixture_ratios = [2, 4, 6, 8, 10]       # O/F

print(f"{'P_cc [Pa]':>12} | {'OF':>4} | {'T_comb [K]':>12}")
print("-"*36)

for P_cc in pressures:
    for OF in mixture_ratios:
        T_cc = C.get_Tcomb(P_cc, OF)
        print(f"{P_cc:12.3e} | {OF:4.1f} | {T_cc:12.3f}")
    print("-"*36)

print(C.get_description())
print(C.get_full_cea_output( Pc=5e5, MR=2, eps=40.0, short_output=1))
