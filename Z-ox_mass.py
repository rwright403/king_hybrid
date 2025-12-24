import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Paths
base_path = "src/results/30_nov_2025_steve"
thrust_path = os.path.join(base_path, "thrust.csv")
m_dot_path = os.path.join(base_path, "m_dot_ox.csv")

# Load CSVs
# Expecting: time(s), value
thrust_data = pd.read_csv(thrust_path)
m_dot_data = pd.read_csv(m_dot_path)

# Extract columns (assumes first column = time, second = value)
t_thrust = thrust_data.iloc[:, 0].to_numpy()
thrust = thrust_data.iloc[:, 1].to_numpy()

# Ox mass flow
t_mdot = m_dot_data.iloc[:, 0].to_numpy()
m_dot_ox = m_dot_data.iloc[:, 1].to_numpy()

# Integrate m_dot to get total mass consumed
m_ox_consumed = np.trapz(m_dot_ox, t_mdot)

print(f"Total oxidizer mass consumed: {m_ox_consumed:.3f} kg")

# Plot both curves
# Plot thrust
plt.figure(figsize=(10, 5))
plt.plot(t_thrust, thrust, label="Thrust (N)")
plt.xlabel("Time (s)")
plt.ylabel("Thrust (N)")
plt.title("Thrust vs Time")
plt.grid(True)
plt.legend()
plt.show()

# Plot m_dot_ox
plt.figure(figsize=(10, 5))
plt.plot(t_mdot, m_dot_ox, label="m_dot_ox (kg/s)")
plt.xlabel("Time (s)")
plt.ylabel("Mass Flow (kg/s)")
plt.title("Ox Mass Flow vs Time")
plt.grid(True)
plt.legend()
plt.show()
