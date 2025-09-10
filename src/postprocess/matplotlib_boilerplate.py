import matplotlib.pyplot as plt

def plot_eq_ox_tank(time, P_ox_tank, m_ox_tank, T_ox_tank):
    plt.subplot(1,3,1)
    plt.plot(time, P_ox_tank, label = "Tank Pressure vs Time")
    plt.xlabel("Time [s]")
    plt.ylabel("Pressure [Pa]")
    plt.title("Ox Tank Pressure")
    plt.legend()
    plt.grid(True)

    plt.subplot(1,3,2)
    plt.scatter(time, m_ox_tank, label = "Tank Mass vs Time")
    plt.xlabel('Time (s)')
    plt.ylabel('Mass (kg)')
    plt.title('Mass vs. Time')
    plt.legend()
    plt.grid(True)

    plt.subplot(1,3,3)
    plt.scatter(time, T_ox_tank, label = "Tank Temp vs Time")
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature vs. Time')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_non_eq_ox_tank(time, P_ox_tank, m_ox_tank_liq, m_ox_tank_gas, T_liq_ox_tank, T_gas_ox_tank, T_sat_ox_tank, T_liq_wall_ox_tank, T_gas_wall_ox_tank):
    plt.subplot(1,3,1)
    plt.scatter(time, P_ox_tank, "model_tank", color = "b" )
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (Pa)')
    plt.title('Pressure vs. Time')
    plt.legend()
    plt.grid(True)


    plt.subplot(1,3,2)
    plt.scatter(time, m_ox_tank_liq, label = "model_liquid", color = "b" )
    plt.scatter(time, m_ox_tank_gas, label = "model_gas", color = "darkorange" )
    plt.xlabel('Time (s)')
    plt.ylabel('Mass (kg)')
    plt.title('Mass vs. Time')
    plt.legend()
    plt.grid(True)


    plt.subplot(1,3,3)
    plt.scatter(time, T_liq_ox_tank, label = "model_liquid", color = "b" )
    plt.scatter(time, T_gas_ox_tank, label = "model_gas", color = "darkorange" )
    plt.scatter(time, T_sat_ox_tank, label = "model_T_sat", color =  "g" )
    plt.scatter(time, T_liq_wall_ox_tank, label = "model_WALL liquid", color = "r" )
    plt.scatter(time, T_gas_wall_ox_tank, label = "model_WALL gas", color = "mediumorchid" )
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature vs. Time')
    plt.legend()
    plt.grid(True)
    plt.show()

#works for both hybrid and liquid
def plot_full_stack(time, P_ox_tank, P_fuel_tank, P_cc, thrust):
        plt.subplot(1,2,1)
        plt.plot(time,thrust)
        plt.xlabel('Time (s)')
        plt.ylabel('Thrust (N)')
        plt.title('Thrust Curve')
        plt.grid(True)

        plt.subplot(1,2,2)
        plt.plot(time,P_ox_tank, label = "P_ox_tank")
        plt.plot(time,P_fuel_tank, label = "P_fuel_tank")
        plt.plot(time,P_cc, label = "P_cc")
        plt.xlabel('Time (s)')
        plt.ylabel('Pressure (Pa)')
        plt.title('System Pressures Over Time')
        plt.grid(True)
        plt.legend()