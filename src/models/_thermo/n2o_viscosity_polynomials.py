### ig polynomials here!!!!
def gas_dynamic_visc_polynomial(T):
    # Polynomial coefficients
    A = 2.1150E-6
    B = 0.46420
    C = 305.70
    D = 0.0

    # Apply temperature limits
    if 182 < T and T < 1000:
        dvisc = A*T**B / ( 1 + C/T + D/T**2)
        return dvisc
    raise ValueError("Temperature outside of function bounds!")

def liq_dynamic_visc_polynomial(T):
    # Polynomial coefficients
    A = 0.001877085
    B = -1.1864E-5  #NOTE: Thesis reported B = 1.1864E-5 but this will return a dynamic viscosity orders of magnitude higher than expected. From testing this function found it was likely signed wrong in the thesis.
    C = 1.928E-8 

    # Apply temperature limits
    if 182 < T and T < 1000:
        dvisc = A + B*T + C*T**2
        return dvisc
    raise ValueError("Temperature outside of function bounds!")