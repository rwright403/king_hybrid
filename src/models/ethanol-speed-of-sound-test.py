import CoolProp.CoolProp as CP

# Example for ethanol
fluid = "Ethanol"
T = 275  # Temperature in K
P = 4e6  # Pressure in Pa

# Calculate speed of sound
speed_of_sound = CP.PropsSI("A", "T", T, "P", P, fluid)
print(f"Speed of sound in ethanol: {speed_of_sound:.2f} m/s")
