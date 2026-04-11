import os
import numpy as np
import matplotlib.pyplot as plt

# =========================
# USER INPUT
# =========================
folder_path = r"src/results/Skvader_config_35_FINAL_T_MAX"  # <-- change this

# File names
p_cc_file = "P_cc.csv"
p_ox_file = "P_ox_tank.csv"
thrust_file = "thrust.csv"

# =========================
# CONSTANTS
# =========================
PA_TO_PSI = 1 / 6894.757293168  # exact Pa → psi conversion

# =========================
# LOAD DATA
# =========================
t_cc, p_cc_pa = np.loadtxt(
    os.path.join(folder_path, p_cc_file),
    delimiter=",",
    skiprows=1,
    unpack=True
)

t_ox, p_ox_pa = np.loadtxt(
    os.path.join(folder_path, p_ox_file),
    delimiter=",",
    skiprows=1,
    unpack=True
)

t_thrust, thrust = np.loadtxt(
    os.path.join(folder_path, thrust_file),
    delimiter=",",
    skiprows=1,
    unpack=True
)

# Convert pressures to psi
p_cc_psi = p_cc_pa * PA_TO_PSI
p_ox_psi = p_ox_pa * PA_TO_PSI

# =========================
# PLOTTING
# =========================
fig, axs = plt.subplots(1, 2, figsize=(12, 5))  # 1 row, 2 columns


# Overall figure title
fig.suptitle("S1 Simulated at T_MAX=31°C", fontsize=16)

# --- Thrust plot ---
axs[0].plot(t_thrust, thrust, color='red')
axs[0].set_xlabel("Time [s]")
axs[0].set_ylabel("Thrust [N]")
axs[0].set_title("Thrust vs Time")
axs[0].grid(True)

# --- Pressure plot ---
axs[1].plot(t_cc, p_cc_psi, label="Chamber Pressure")
axs[1].plot(t_ox, p_ox_psi, label="Ox Tank Pressure")
axs[1].set_xlabel("Time [s]")
axs[1].set_ylabel("Pressure [Psi]")
axs[1].set_title("Pressures vs Time")
axs[1].legend()
axs[1].grid(True)

plt.tight_layout()
plt.show()
