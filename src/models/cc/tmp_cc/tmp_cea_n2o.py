from rocketcea.cea_obj_w_units import CEA_Obj, add_new_fuel

# Define nitrous oxide as a "fuel" or "oxidizer"
fuel_properties = """#fuel paraffin  N 2   O 1    wt%=100.00
"""
add_new_fuel('N2O', fuel_properties)

# Create a CEA object with nitrous oxide as both fuel and oxidizer,
# but we’ll only use it in “oxidizer mode”
n2oCEA = CEA_Obj(oxName='N2O', fuelName='N2O')

# Now we can query thermo properties
T = 298.15
P = 30.0  # bar
eps = 1.0  # not used here, but required by API

# Trick: use a giant O/F ratio so it's "all oxidizer"
OF = 1

h = n2oCEA.get_Enthalpies(T, P, OF, eps)[0]

print("h_test: ", h, "at (T, P): ", T, P)
