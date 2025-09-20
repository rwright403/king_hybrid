def system_rhs(t, y, tank, chamber):
    # unpack
    y_tank = y[tank_slice]
    y_cc   = y[cc_slice]

    # tank ODE
    dydt_tank, outputs_tank = tank.ode_system(t, y_tank)

    # chamber ODE (using tank outputs as inputs)
    m_dot_ox = outputs_tank["m_dot_ox"]
    dydt_cc, outputs_cc = chamber.cc_ode_system(t, y_cc, m_dot_ox)

    # combine
    dydt = np.zeros_like(y)
    dydt[tank_slice] = dydt_tank
    dydt[cc_slice]   = dydt_cc

    return dydt

sol = solve_ivp(
    lambda t,y: system_rhs(t,y,tank,chamber),
    t_span=(0, sim_time),
    y0=full_y0,
    method="RK45",
    rtol=1e-6, atol=1e-9
)
