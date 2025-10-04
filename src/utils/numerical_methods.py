import numpy as np

"""
def secant(func, x1):
    x_eps = x1 * 0.005  # Set the tolerance to be 0.5% of init guess
    x2 = x1 -x1 * 0.01  # Set a second point 1% away from the original guess
    F1 = func(x1)  # Evaluate function at x1
    F2 = func(x2)  # Evaluate function at x2
    kk = 1  # Set up counter
    kk_max = 1000

    while np.abs(x2 - x1) >= (x_eps) and kk < kk_max:  # While error is too large and counter is less than max
        x3 = x2 - (F2 * (x2 - x1) / (F2 - F1)) 
        x1 = x2  # Move everything forward
        x2 = x3
        F1 = F2
        F2 = func(x2) 
        if (F1 == F2):
            return x2
        kk = kk + 1
    x = x2
    return x
"""
"""
def rk4_step(ode_func, t, y, dt, *args):
    k1 = ode_func(t, y, *args)
    y2 = [y_i + 0.5*dt*k1_i for y_i, k1_i in zip(y, k1)]
    k2 = ode_func(t + 0.5*dt, y2, *args)

    y3 = [y_i + 0.5*dt*k2_i for y_i, k2_i in zip(y, k2)]
    k3 = ode_func(t + 0.5*dt, y3, *args)

    y4 = [y_i + dt*k3_i for y_i, k3_i in zip(y, k3)]
    k4 = ode_func(t + dt, y4, *args)

    y_new = [y_i + (dt/6.0)*(k1_i + 2*k2_i + 2*k3_i + k4_i)
             for y_i, k1_i, k2_i, k3_i, k4_i in zip(y, k1, k2, k3, k4)]

    return y_new
"""

"""
This exists because I made bad decisions that seemed ok at the time:
My development of this script involved testing a lot of models, so I decided to make
them loosely coupled so I could quickly try things.
Because they are loosely coupled, I couldn't use any RK implementation with adaptive timestepping
I only ended up using really 2 different models. My original algo was really slow so I was planning
to use the equilibrium tank more, but I refactored my algo and vectorized some things to improve speed.
This made the previous idea redundant.
"""
def rk4_step(ode_func, t, y, dt, *args):
    y = np.asarray(y, dtype=float)

    k1 = np.asarray(ode_func(t, y, *args))
    k2 = np.asarray(ode_func(t + 0.5*dt, y + 0.5*dt*k1, *args))
    k3 = np.asarray(ode_func(t + 0.5*dt, y + 0.5*dt*k2, *args))
    k4 = np.asarray(ode_func(t + dt, y + dt*k3, *args))

    return y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
