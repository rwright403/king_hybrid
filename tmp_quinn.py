import CoolProp.CoolProp as CP
import math

pi = math.pi

#----------------------------
#N2O Liquid Density Function
#----------------------------
def n2o_liquid_density(T):

    rho = CP.PropsSI("D", "T", T, "Q", 0, "NitrousOxide")

    return rho

def n2o_vapor_density(T):
    return CP.PropsSI("D", "T", T, "Q", 1, "NitrousOxide")

#----------------------------
#Hollow Cylinder Mass Function
#----------------------------
def hollow_cylinder_mass(rho, L, d_inner, d_outer):
    return rho * pi * L * ((d_outer / 2)**2 - (d_inner / 2)**2)

#----------------------------
#Tank and Oxidizer Calculation
#----------------------------
#Given values
v_ox_tank = 0.029 #* 0.85       # m^3
rho_al = 2700                  # kg/m^3
d_outer = 0.1413                # m
d_inner = 0.1223                # m

#Ox tank length from volume & inner radius
r_inner = d_inner / 2
L = v_ox_tank / (pi * r_inner**2)

#Compute masses
cylinder_m = hollow_cylinder_mass(rho_al, L, d_inner, d_outer)

#N2O density
T = 273.15 + 31#17   # 18°C
#P = 4.7e6         # Pa
rho_ox = n2o_liquid_density(T)
rho_vap = n2o_vapor_density(T)

actual_rho = 17.7/v_ox_tank
print("actual rho: ", actual_rho,)


#nitrous_m = v_ox_tank * rho_ox * 0.9 # 90% fill
nitrous_m = (0.85*v_ox_tank)*rho_ox + (0.15*v_ox_tank)*rho_vap # 90% fill

print("sol rho: ", nitrous_m/v_ox_tank)
wetted_tank = cylinder_m + nitrous_m

#----------------------------
#Combustion Chamber (Hybrid Fuel Grain)
#----------------------------
OF = 7.5
fuel_m = nitrous_m / OF
rho_fuel = 900.0
V_fuel = fuel_m / rho_fuel
A_inner = pi *r_inner *r_inner
A_port = 0.0029
L_port = V_fuel / (A_inner - A_port)

#----------------------------
#Print Results
#----------------------------
print(f"Ox_Tank_Length: {L:.4f} m")
print(f"Ox_Tank_Mass: {cylinder_m:.2f} kg")
print(f"Nitrous_Mass: {nitrous_m:.2f} kg")
print(f"Wetted_Tank_Mass: {wetted_tank:.2f} kg")
print(f"N2O Liquid Density: {rho_ox:.2f} kg/m^3")
print(f"Fuel_Mass: {fuel_m:.2f} kg")
print(f"Combustion_Chamber_Length: {L_port:.4f} m")