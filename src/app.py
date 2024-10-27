import argparse
import importlib
from src.thrust_curve import thrust_curve

#from src.sensitivity_analysis import sensitivity_analysis
from src.model_validation import model_validation
from src.model_validation import nitrous_validation
from src.prelim_wizard import prelim_wizard
from src.sensitivity_analysis import sensitivity_analysis



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
        print("6 --> Oxidizer Tank and Injector Validation")

        user_input = input("Enter number to select analysis: ")

    if user_input =='1':
        #print("Running Thrust Curve from inputs in constants.py") #TODO: UPDATE
        thrust_curve.run_thrust_curve(program_input)

    if user_input =='2':
        #print("Running Flight Sim from inputs in constants.py") #TODO: UPDATE
        thrust_curve.run_thrust_curve(program_input)
        from src.flight_sim import flight_sim

    if user_input =='3':
        sensitivity_analysis.run_sensitivity_analysis(program_input)

    if user_input =='4':
        thrust_curve.run_thrust_curve(program_input)
        model_validation.validate(program_input)

    if user_input == '5':
        prelim_wizard.magic(program_input)
    
    if user_input == '6':
        nitrous_validation.run_nitrous_validation(program_input)