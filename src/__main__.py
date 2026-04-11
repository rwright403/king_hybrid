import os
import argparse
import importlib
import time

from src.utils.kwargs_builder import build_kwargs
from src.utils.flight_sim_kwargs_builder import build_flight_sim_kwargs
from src.sim.prop_sim import prop_sim
from src.sim.flight_sim import flight_sim
from src.models.aerostruct.aerostruct import aerostruct
from src.prelim_wizard import prelim_wizard
from src.sensitivity_analysis import sensitivity_analysis
from src.postprocess.engine_overview import engine_overview
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

        start_time = time.time()
        prop_results = prop_sim(kwargs)
        end_time = time.time()
        elapsed = end_time - start_time

        print(f"\nSimulation completed in real time {elapsed:.2f} seconds ({elapsed/60:.2f} min)\n")

        # Make a directory path for this case
        case_dir = os.path.join("src/results", str(input_file))

        # Save outputs into that folder
        build_rocketpy_input_csv(prop_results, "m_dot_ox", output_dir=case_dir)
        build_rocketpy_input_csv(prop_results, "P_ox_tank", output_dir=case_dir)
        
        #fuel tank run
        if kwargs["models_kwargs"]["fuel_tank_model"] is not None:      # liquid engine
            build_rocketpy_input_csv(prop_results, "m_dot_fuel", output_dir=case_dir)
        if kwargs["models_kwargs"]["cc_model"] != 3:
            build_rocketpy_input_csv(prop_results, "thrust", output_dir=case_dir)
            build_rocketpy_input_csv(prop_results, "m_dot_cc", output_dir=case_dir)
        
        # hybrid engine
        if kwargs["models_kwargs"]["cc_model"] == 1:           
            build_rocketpy_input_csv(prop_results, "OF", output_dir=case_dir)

            # tmp add nozzle stuff here, todo redistribute to cc
            build_rocketpy_input_csv(prop_results, "P_cc", output_dir=case_dir)
            build_rocketpy_input_csv(prop_results, "v_exit", output_dir=case_dir)
            build_rocketpy_input_csv(prop_results, "P_exit", output_dir=case_dir)


        engine_overview(program_input, prop_results, program_input.mode)
        plot_sim_results(program_input, prop_results, program_input.mode, program_input.save_path)

    if user_input =='2':
        flight_sim_kwargs, mass_data = build_flight_sim_kwargs(input_file, program_input)
        rocket, flight = flight_sim(flight_sim_kwargs)

        
        if flight_sim_kwargs["rocketpy_models_kwargs"]["mass_model"] == 2:           # this mass model gives us mass distribution, can sol shear + bending
            rkt_len = flight_sim_kwargs["rocketpy_rocket_kwargs"]["nose_position"]
            #aerostruct(rkt_len, mass_data, rocket, flight)

            #just add nozzle study here! - tmp, TODO: MOVE THIS!!!!
            from src.postprocess.nozzle_performance import nozzle_performance
            import numpy as np
            
            
            #write rocketpy flight atmospheric pressures to csv:
            import pandas as pd
            
            t = flight.time

            # Evaluate pressure at each time
            P_atm = [flight.env.pressure(ti) for ti in t]
            Alt = [flight.z(ti) for ti in t]

            # Make a DataFrame
            df = pd.DataFrame({
                "time": flight.time,               # seconds
                "P_atm": P_atm,  # Pa or kPa depending on RocketPy units
                "Alt": Alt,
            })

            # Write CSV
            df.to_csv("atm_pressure.csv", index=False)
            case_dir = os.path.join("src/results", str(input_file))
            build_rocketpy_input_csv(df, "P_atm", output_dir=case_dir)

            df.to_csv("alt.csv", index=False)
            case_dir = os.path.join("src/results", str(input_file))
            build_rocketpy_input_csv(df, "Alt", output_dir=case_dir)

            df.to_csv("P_ox_tank.csv", index=False)
            case_dir = os.path.join("src/results", str(input_file))
            build_rocketpy_input_csv(df, "P_ox_tank", output_dir=case_dir)


            A_throat = 0.25*np.pi*(program_input.d_throat**2)
            A_exit = A_throat*program_input.expratio

            nozzle_performance(
                ve_csv=os.path.join("src/results", str(input_file), "v_exit.csv"),
                mdot_csv=os.path.join("src/results", str(input_file), "m_dot_cc.csv"),
                pc_csv=os.path.join("src/results", str(input_file), "P_cc.csv"),
                pe_csv=os.path.join("src/results", str(input_file), "P_exit.csv"),
                pa_csv=os.path.join("src/results", str(input_file), "P_atm.csv"),
                alt_csv=os.path.join("src/results", str(input_file), "Alt.csv"),
                At=A_throat,
                Ae=A_exit,
                t_burn=0.8,
            )


        

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
