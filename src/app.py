import argparse
import importlib
from src.sim.sim import sim
from src.prelim_wizard import prelim_wizard
from src.sensitivity_analysis import sensitivity_analysis
from src.postprocess.plot_sim_results import plot_sim_results



def run(input_file):

    # Dynamically import the input file
    program_input = importlib.import_module(f"src.inputs.{input_file}")


    user_input = 0
    while(user_input ==0):
        print("\n")
        print("1 --> Sim")
        print("2 --> Flight Sim")
        """print("3 --> Sensitivity Analysis")
        print("4 --> Model Validation")
        print("5 --> Prelim Design Wizard")
        print("6 --> Oxidizer Tank and Injector Validation")"""

        user_input = input("Enter number to select analysis: ")

    if user_input =='1':
        results = sim(program_input)
        plot_sim_results(program_input, results, program_input.mode, program_input.save_path)

    if user_input =='2':
        sim.sim(program_input)
#NOTE: THIS ONLY WORKS IF THE SIM IS A FULL STACK THRUST CURVE
        from src.flight_sim import flight_sim

    if user_input =='3':
        sensitivity_analysis.run_sensitivity_analysis(program_input)

    if user_input == '5':
        prelim_wizard.magic(program_input)