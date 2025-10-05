import argparse
import importlib
from src.utils.kwargs_builder import build_kwargs
from src.utils.flight_sim_kwargs_builder import build_flight_sim_kwargs
from src.sim.prop_sim import prop_sim
from src.sim.flight_sim import flight_sim
from src.prelim_wizard import prelim_wizard
from src.sensitivity_analysis import sensitivity_analysis
from src.postprocess.plot_sim_results import plot_sim_results
from src.utils.build_rocketpy_input_csv import build_rocketpy_input_csv


def run(input_file):

    """
    print("WARNING, THE FOLLOWING TEST CODE HAS BEEN ADDED FOR DEBUGGING:\n"
    "")
    """

    # Dynamically import the input file
    program_input = importlib.import_module(f"src.inputs.{input_file}")


    user_input = 0
    while(user_input ==0):
        print("\n")
        print("1 --> Engine Sim")
        print("2 --> Flight Sim")
        #("3 --> Sensitivity Analysis") #TODO: refactor
        print("5 --> Prelim Design Wizard")

        user_input = input("Enter number to select analysis: ")

    if user_input =='1':

        kwargs = build_kwargs(program_input)
        prop_results = prop_sim(kwargs)
        
        build_rocketpy_input_csv(prop_results, "m_dot_ox")
        build_rocketpy_input_csv(prop_results, "m_dot_fuel")
        build_rocketpy_input_csv(prop_results, "thrust")

        plot_sim_results(program_input, prop_results, program_input.mode, program_input.save_path)

    if user_input =='2':
        flight_sim_kwargs = build_flight_sim_kwargs(program_input)
        flight = flight_sim(flight_sim_kwargs)

    """
    #TODO: refactor
    if user_input =='3':
        sensitivity_analysis.run_sensitivity_analysis(program_input)
    """

    if user_input == '5':
        prelim_wizard.magic(program_input)

        print(program_input)



if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Run the simulation with a specified input file.')
    parser.add_argument('input_file', type=str, help='The name of the input file to use (e.g., constants1) without the .py extension')
    
    args = parser.parse_args()
    
    # Run the program with the specified input file
    run(args.input_file)