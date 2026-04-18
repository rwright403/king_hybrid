import numpy as np
from scipy.optimize import root


def residual(x_new, y_old, dt, k):
    """
    x_new = [y_new, p_new]
    """
    y_new, p_new = x_new

    R1 = y_new - y_old - dt * (-k * y_new + p_new) #deriv
    R2 = p_new - y_new**2 #algebra

    return np.array([R1, R2])


# settings
dt = 0.1
k = 2.0
n_steps = 20

# initial condition
y = 1.0
p = y**2

history = [(0.0, y, p)]

for n in range(n_steps):
    # good initial guess = previous converged state
    x_guess = np.array([y, p])

    sol = root(residual, x_guess, args=(y, dt, k), method="hybr")

    if not sol.success:
        raise RuntimeError(f"Nonlinear solve failed at step {n}: {sol.message}")

    y, p = sol.x
    t = (n + 1) * dt
    history.append((t, y, p))

# print results
for t, y, p in history:
    print(f"t = {t:5.2f}   y = {y: .6f}   p = {p: .6f}")