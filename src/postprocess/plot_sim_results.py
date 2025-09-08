from src.postprocess.matplotlib_boilerplate import *

#TODO: if validation data exists and mode is validation data, then u should graph it 
def plot_sim_results(inputs, results: dict, mode: str, save_path: str = None):

    time = results["time"]



    if mode == "ox_tank":

        if inputs.ox_tank_model == 1: #equilibrium tank, plot P_tank, T_tank, m_tank?
            plot_eq_ox_tank(time, results["P_ox_tank"], results["m_ox_tank"], results["T_ox_tank"])
            
        if inputs.ox_tank_model == 2: #non-equilibrium tank, plot P_tank, P_sat, T_liq, T_gas, T_sat, 
            plot_non_eq_ox_tank(time, results["P_ox_tank"], results["m_ox_tank_liq"], results["m_ox_tank_gas"], results["T_liq_ox_tank"], results["T_gas_ox_tank"], results["T_sat_ox_tank"], results["T_liq_wall_ox_tank"], results["T_gas_wall_ox_tank"])
            


    elif mode == "fuel_tank":
        plt.plot(time, results["P_fuel_tank"], label="P_fuel Tank")
        plt.xlabel("Time [s]")
        plt.ylabel("Pressure [Pa]")
        plt.title("Fuel Tank Pressure")
        plt.legend()
        plt.grid(True)

    elif mode == "full_stack":
        if len(results["P_fuel_tank"]) > 0: # P_tank added if this col has results
            plot_liquid_full_stack(time, results["P_ox_tank"], results["P_fuel_tank"], results["P_cc"], results["thrust"])
        else:
            plot_hybrid_full_stack(time, results["P_ox_tank"], results["P_cc"], results["thrust"])


    # Save/show logic
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    plt.show()
