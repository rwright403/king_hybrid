from rocketcea.cea_obj import CEA_Obj

def test_propellant(name):
    try:
        # Just try to create a CEA object with the propellant
        cea = CEA_Obj(oxName="LOX", fuelName=name)
        # Query a simple property to force CEA to check validity
        isp = cea.get_Isp(Pc=1000, MR=2.5, eps=10)
        print(f"✅ '{name}' is valid. Example Isp: {isp:.2f} s")
    except Exception as e:
        print(f"❌ '{name}' is NOT valid. Error: {e}")

if __name__ == "__main__":
    test_propellant("Ethanol")
