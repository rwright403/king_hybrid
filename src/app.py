import argparse
import importlib
from src.thrust_curve import thrust_curve
#from src.sensitivity_analysis import sensitivity_analysis



def run(input_file):

    # Dynamically import the input file
    program_input = importlib.import_module(f"src.inputs.{input_file}")


    user_input = 0
    while(user_input ==0):
        print("\n")
        print("1 --> Thrust Curve")
        print("2 --> Flight Sim")
        print("3 --> Sensitivity Analysis")
        print("4 --> Model Validation")
        print("5 --> Prelim Design Wizard")

        user_input = input("Enter number to select analysis: ")

    if user_input =='1':
        print("Running Thrust Curve from inputs in constants.py")
        thrust_curve.run_thrust_curve(program_input)

    if user_input =='2':
        print("Running Flight Sim from inputs in constants.py")
        thrust_curve.run_thrust_curve(program_input)
        from src.flight_sim import flight_sim

    if user_input =='3':
        pass
        #sensitivity_analysis.sensitivityAnalysis(program_input)

    if user_input =='4':
        thrust_curve.run_thrust_curve(program_input)
        from src.model_validation import model_validation

    if user_input == '5':
        from src.prelim_wizard import prelim_wizard