from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
import csv

# Recovery Parameters

mn          = 3.45                  # Mass of Upper Airframe                [kg]
mr          = 10.94                 # Mass of Lower Airframe                [kg]
m_parachute = 0.07                  # Mass of Parachute                     [kg]
k           = 4942                  # Shock Cord Spring Constant            [N/m]
cd          = 1.3                   # Coefficient of drag                   [-]
proj_area   = 0.3                   # Projected Parachute Area              [m^2]
n_lines     = 8                     # Number of Shroud Lines                [-] 
line_angle  = np.deg2rad(30)        # Shrould Line Angle                    [Deg]

# Environmental Parameters

rho         = 0.38                  # Air Density at Deployment             [kg/m^3]
g           = 9.81                  # Gravity                               [m/s^2]
v_wind      = 55.6                  # Wind Velocity at Deployment           [m/s]

# Dynamic Parameters

xp          = 0                     # Parachute Initial Position            [m]     {always zero}
xn          = 0                     # Upper Airframe Initial Position       [m]     {always zero}
xr          = 0                     # Lower Airframe Initial Position       [m]     {always zero}

xp_dot      = 0                     # Parachute Initial Velocity            [m/s]   {zero unless reefed chute}
xn_dot      = 40.23                 # Initial Velocity of Upper Airframe    [m/s]
xr_dot      = 40.23                 # Initial Velocity of Lower Airframe    [m/s]

# Time Vector

t_start     = 0                     # Lower Time Duration to Solve Over     [s]
t_end       = 10                    # Upper Time Duration to Solve Over     [s]
nt          = 1000                  # Number of Intervals                   [-]     {larger value = more refined data}

t = np.linspace(t_start, t_end, nt) # Time Vector

# Initial Conditions

y0 = [xp, xn, xr, xp_dot, xn_dot, xr_dot]

################################################################################
############################### Computation ####################################
################################################################################

# Calculated Variables

c = 0.5 * cd * rho * proj_area

# Differential Function
def dy_dt(t, y, m_parachute, mn, mr, c, k, g, v_wind):
    
    xp, xn, xr, xp_dot, xn_dot, xr_dot = y

    # First-order equations
    dxp_dt = xp_dot
    dxn_dt = xn_dot
    dxr_dt = xr_dot

    # Second-order equations
    dxp_dotdt = (m_parachute*g - c*(xp_dot + v_wind)**2 - k*(xp - xr) - k*(xp - xn)) / m_parachute
    dxn_dotdt = (mn*g - k*(xn - xp)) / mn
    dxr_dotdt = (mr*g - k*(xr - xp)) / mr

    return [dxp_dt, dxn_dt, dxr_dt, dxp_dotdt, dxn_dotdt, dxr_dotdt]

# Solve 3 Second Order ODEs
solution = solve_ivp(dy_dt, (t_start, t_end), y0, args=(m_parachute, mn, mr, c, k, g, v_wind), t_eval=t)

# ---------------------- Prepare data for plotting ----------------------
t_sol = solution.t
xp     = solution.y[0]
xn     = solution.y[1]
xr     = solution.y[2]
xp_dot = solution.y[3]
xn_dot = solution.y[4]
xr_dot = solution.y[5]

# Shock cord tensions
upper_shock_tension = k*(xn - xp)
lower_shock_tension = k*(xr - xp)
shock_cord_total = upper_shock_tension + lower_shock_tension

# Parachute acceleration (with adaptive time step)
parachute_acceleration = np.diff(xp_dot) / np.diff(t_sol)
F_vertical = (m_parachute) * parachute_acceleration + c * (v_wind + xp_dot[:-1])**2

# Shroud line tensions
shroud_line_tensions = np.zeros((len(F_vertical), n_lines))
for i in range(n_lines):
    shroud_line_tensions[:, i] = F_vertical / (n_lines * np.cos(line_angle))

xp_accel = np.diff(xp_dot) / np.diff(t_sol)
xn_accel = np.diff(xn_dot) / np.diff(t_sol)
xr_accel = np.diff(xr_dot) / np.diff(t_sol)

# ---------------------- Acceleration Graphs ----------------------
plt.figure(figsize=(10,6))
plt.plot(t_sol[:-1], xp_accel, label='Parachute Acceleration')
plt.plot(t_sol[:-1], xn_accel, label='Upper Airframe Acceleration')
plt.plot(t_sol[:-1], xr_accel, label='Lower Airframe Acceleration')
plt.xlabel('Time [s]')
plt.ylabel('Acceleration [m/s²]')
plt.title('Accelerations vs Time')
plt.legend()
plt.grid(True)
plt.savefig('acceleration-vs-time.png')
plt.show()

# ---------------------- Position plot ----------------------
plt.figure(figsize=(10,6))
plt.plot(t_sol, xp, label='Parachute Position')
plt.plot(t_sol, xn, label='Upper Airframe Position')
plt.plot(t_sol, xr, label='Lower Airframe Position')
plt.xlabel('Time [s]')
plt.ylabel('Displacement [m]')
plt.title('Parachute and Airframe Positions vs Time')
plt.legend()
plt.grid(True)
plt.savefig("disp-vs-t.png")
plt.show()

# ---------------------- Shock cord tensions ----------------------
plt.figure(figsize=(10,6))
plt.plot(t_sol, upper_shock_tension, label='Upper Shock Cord')
plt.plot(t_sol, lower_shock_tension, label='Lower Shock Cord')
plt.plot(t_sol, shock_cord_total, '--', label='Total Shock Cord')
plt.xlabel('Time [s]')
plt.ylabel('Force [N]')
plt.title('Shock Cord Forces vs Time')
plt.legend()
plt.grid(True)
plt.savefig('shockcord-force-vs-time.png')
plt.show()

################################################################################
################################# RESULTS ######################################
################################################################################

print("\nUsing the following solver parameters:\n")

print(f"Upper Airframe Mass:\t\t {mn:.2f} [kg]")
print(f"Lower Airframe Mass:\t\t {mr:.2f} [kg]")
print(f"Parachute Mass:\t\t {m_parachute:.2f} [kg]")
print(f"Shock Cord k:\t\t {k:.0f} [N/m]")
print(f"Air Density:\t\t {rho:.2f} [kg/m^3]")
print(f"Wind Speed:\t\t {v_wind:.1f} [m/s]")
print(f"Drag Coefficient:\t {cd:.2f}")
print(f"Projected Parachute Area:\t\t {proj_area:.2f} [m^2]")
print(f"Upper Airframe Velocity at inflation:\t {xn_dot[0]:.2f} [m/s]")
print(f"Lower Airframe Velocity at Inflation:\t {xr_dot[0]:.2f} [m/s]")
print(f"Parachute velocity at inflation:\t {xp_dot[0]:.2f} [m/s]")

print("\nProduced the following results:\n")

print(f"Parachute Damping Coeff:\t\t {c:.2f} [N/m/s]")
print(f"Maximum Upper Shock Cord Tension:\t {max(upper_shock_tension):.0f} [N]")
print(f"Maximum Lower Shock Cord Tension:\t {max(lower_shock_tension):.0f} [N]")
print(f"Maximum Shroud Line Tension:\t\t {np.max(shroud_line_tensions):.0f} [N]")
print(f"Maximum Upper Airframe G load:\t\t {max(upper_shock_tension)/(mn*g):.2f} [G]")
print(f"Maximum Lower Airframe G load:\t\t {max(lower_shock_tension)/(mr*g):.2f} [G]")