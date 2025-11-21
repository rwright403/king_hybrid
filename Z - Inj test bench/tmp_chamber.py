import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from mpl_toolkits.mplot3d import Axes3D

# --- Parameters ---
fluid = 'NitrousOxide'   # Or 'N2O' if you have aliases set up
V_chamber = 0.000510994        # m^3, example downstream chamber volume
T_range = np.linspace(260, 300, 30)  # K, range of temperatures
m_range = np.linspace(0.001, 0.1, 30)  # kg, mass of fluid in chamber

# --- Arrays for plotting ---
T_grid, m_grid = np.meshgrid(T_range, m_range)
P_grid = np.zeros_like(T_grid)

# --- Compute pressure for each T, m ---
for i in range(T_grid.shape[0]):
    for j in range(T_grid.shape[1]):
        T = T_grid[i, j]
        m = m_grid[i, j]
        rho = m / V_chamber
        try:
            P = PropsSI('P', 'D', rho, 'T', T, fluid)
        except:
            P = np.nan
        P_grid[i, j] = P / 1e5  # convert to bar

# --- 3D Plot ---
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(m_grid*1000, T_grid, P_grid, cmap='viridis', alpha=0.9)

# --- Add a horizontal plane at 30 bar ---
P_target = 30  # bar
ax.plot_surface(
    m_grid*1000,  # x (mass)
    T_grid,       # y (temperature)
    np.full_like(P_grid, P_target),  # z constant at 30 bar
    color='red',
    alpha=0.3,
    label='30 bar plane'
)

# --- Labels ---
ax.set_xlabel('Mass [g]')
ax.set_ylabel('Temperature [K]')
ax.set_zlabel('Pressure [bar]')
ax.set_title(f'{fluid} Pressure vs. Mass & Temperature\nV_chamber = {V_chamber*1e3:.1f} L')

plt.tight_layout()
plt.show()