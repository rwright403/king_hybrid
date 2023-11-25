""" HYBRID 2024!!!! """



def run():
    user_input = 0
    while(user_input ==0):
        print("\n")
        print("1 --> Thrust Curve")
        print("2 --> Flight Sim")
        print("3 --> Sensitivity Analysis")
        print("4 --> Model Validation")

        user_input = input("Enter number to select analysis:")

    if user_input =='1':
        print("Running Thrust Curve from inputs in constants.py")
        from src.thrust_curve import thrust_curve

    if user_input =='2':
        print("Running Flight Sim from inputs in constants.py")
        from src.thrust_curve import thrust_curve
        from src.flight_sim import flight_sim

    if user_input =='3':
        exit()

    if user_input =='4':
        from src.thrust_curve import thrust_curve
        from src.model_validation import model_validation