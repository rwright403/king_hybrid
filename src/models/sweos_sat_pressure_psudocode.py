import numpy as np
from scipy.optimize import brentq

import matplotlib.pyplot as plt

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

def P_sat_anc(T):
    # Polynomial coefficients
    A = 4.80716087
    B = 967.819748
    C = 19.6368887

    if 140 < T and T < 310:
        P_sat_est = 100000*( 10**(A-(B/(T+C))) )
        return P_sat_est
    raise ValueError("Temperature outside of function bounds!")

def T_sat_anc(P):
    # Polynomial coefficients
    A = 4.80716087
    B = 967.819748
    C = 19.6368887

    T_sat_est = B/(A - np.log10(P/100000)) - C

    if T_sat_est > 140 and 310 < T_sat_est:
        raise ValueError("Temperature outside of function bounds!")
    return T_sat_est


def explicit_helmholtz_derivs(rho, T):

    # Constants for N2O
    R = 8.3144598 / 44.0128 * 1000 # Specific Gas constant (J/kg*K)
    T_c = 309.52  # Critical Temperature (K)
    rho_c = 452.0115  # Critical Density (kg/m^3)

    n0 = np.array([0.88045, -2.4235, 0.38237, 0.068917, 0.00020367, 0.13122, 0.46032,
        -0.0036985, -0.23263, -0.00042859, -0.042810, -0.023038])
    n1 = n0[0:5]
    n2 = n0[5:12]
    a1 = 10.7927224829
    a2 = -8.2418318753
    c0 = 3.5
    v0 = np.array([2.1769, 1.6145, 0.48393])
    u0 = np.array([879, 2372, 5447])
    t0 = np.array([0.25, 1.125, 1.5, 0.25, 0.875, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5])
    d0 = np.array([1, 1, 1, 3, 7, 1, 2, 5, 1, 1, 4, 2])
    P0 = np.array([1, 1, 1, 2, 2, 2, 3])
    t1 = t0[0:5]
    t2 = t0[5:12]
    d1 = d0[0:5]
    d2 = d0[5:12]

    # Calculate non-dimensional variables
    tau = T_c / T
    delta = rho / rho_c

    # Calculate explicit Helmholtz energy and derivatives
    ao = a1 + a2 * tau + np.log(delta) + (c0 - 1) * np.log(tau) + np.sum(v0 * np.log(1 - np.exp(-u0 * tau / T_c)))
    ar = np.sum(n1 * tau**t1 * delta**d1) + np.sum(n2 * tau**t2 * delta**d2 * np.exp(-delta**P0))
    ao_tau = a2 + (c0 - 1) / tau + np.sum(v0 * u0 / T_c * np.exp(-u0 * tau / T_c) / (1 - np.exp(-u0 * tau / T_c)))
    ao_tautau = -(c0 - 1) / tau**2 + np.sum(-v0 * u0**2 / T_c**2 * np.exp(-u0 * tau / T_c) / (1 - np.exp(-u0 * tau / T_c))**2)
    ar_tau = np.sum(n1 * t1 * tau**(t1 - 1) * delta**d1) + np.sum(n2 * t2 * tau**(t2 - 1) * delta**d2 * np.exp(-delta**P0))
    ar_tautau = np.sum(n1 * t1 * (t1 - 1) * tau**(t1 - 2) * delta**d1) + np.sum(n2 * t2 * (t2 - 2) * tau**(t2 - 2) * delta**d2 * np.exp(-delta**P0))
    ar_delta = np.sum(n1 * d1 * delta**(d1 - 1) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))
    ar_deltadelta = np.sum(n1 * d1 * (d1 - 1) * delta**(d1 - 2) * tau**t1) + np.sum(n2 * tau**t2 * delta**(d2 - 2) * ((d2 - P0 * delta**P0) * (d2 - 1 - P0 * delta**P0) - P0**2 * delta**P0) * np.exp(-delta**P0))
    ar_deltatau = np.sum(n1 * d1 * t1 * delta**(d1 - 1) * tau**(t1 - 1)) + np.sum(n2 * t2 * tau**(t2 - 1) * delta**(d2 - 1) * (d2 - P0 * delta**P0) * np.exp(-delta**P0))

    return {
        'R': R,
        'tau': tau,
        'delta': delta,
        'ao': ao,
        'ar': ar,
        'ao_tau': ao_tau,
        'ar_tau': ar_tau,
        'ao_tautau': ao_tautau,
        'ar_tautau': ar_tautau,
        'ar_delta': ar_delta,
        'ar_deltadelta': ar_deltadelta,
        'ar_deltatau': ar_deltatau,
    }


def verbose_find_all_sweos_density_roots(T, P_target, rho_min, rho_max, n_points=1000):
    
    rhos = np.linspace(rho_min, rho_max, n_points)
    residuals = [thermo_span_wagner(rho, T, 'p') - P_target for rho in rhos]
    
    signs = np.sign(residuals)
    zero_crossings = np.where(np.diff(signs))[0]
    
     # 🔍 Plot pressure residuals
    plt.figure(figsize=(8, 5))
    plt.plot(rhos, residuals, label='P_calc - P_target')
    plt.axhline(0, color='red', linestyle='--', label='Zero residual')
    plt.xlabel('Density [kg/m³]')
    plt.ylabel('Pressure Residual [Pa]')
    plt.title(f'Residuals of Helmholtz EOS Pressure at T = {T:.2f} K')
    plt.grid(True)
    plt.legend()
    plt.show()

    roots = []
    for idx in zero_crossings:
        try:
            root = brentq(lambda r: thermo_span_wagner(rho=r, T=T, param='p') - P_target,
                          rhos[idx], rhos[idx+1])
            roots.append(root)
        except ValueError:
            pass  # in case brentq fails due to no root in interval
    return roots


def find_all_sweos_density_roots(T, P_target, rho_min, rho_max, n_points):
    
    rhos = np.linspace(rho_min, rho_max, n_points)
    residuals = [thermo_span_wagner(rho=rho, T=T, param='p') - P_target for rho in rhos]
    
    signs = np.sign(residuals)
    zero_crossings = np.where(np.diff(signs))[0]

    roots = []
    for idx in zero_crossings:
        try:
            root = brentq(lambda r: thermo_span_wagner(rho=r, T=T, param='p') - P_target,
                          rhos[idx], rhos[idx+1])
            roots.append(root)
        except ValueError:
            pass  # in case brentq fails due to no root in interval

    if len(roots) < 2:
        raise ValueError("Less than two density roots found")
    return sorted(roots) #returns in order of smallest to largest. NOTE: rho_gas = roots[0] // rho_liq = roots[-1] 

def span2000_residual(P, T):
    # Solve for liquid and vapor densities at this P, T
    roots = find_all_sweos_density_roots(T, P, rho_min=0.1, rho_max=1200,n_points=1000)
    
    rho_gas = roots[0]      # Lowest density (vapor)
    rho_liq = roots[-1]     # Highest density (liquid)

    # Compute residual combining pressure and chemical potential constraints
    prop_gas = explicit_helmholtz_derivs(rho_gas, T)
    prop_liq = explicit_helmholtz_derivs(rho_liq, T)
    RES = P/(prop_liq['R']*T) * (1/rho_gas-1/rho_liq) - np.log(rho_liq/rho_gas) - (prop_liq['ar']-prop_gas['ar'])

    return RES


def sol_init_guess(est, pcnt):
    est_min = (1-pcnt)*est
    est_max = (1+pcnt)*est
    return est_min, est_max

def thermo_span_wagner(rho=None, T=None, P=None, param=None):
    
    if param == 'p':  # Pressure (Pa)
        prop = explicit_helmholtz_derivs(rho, T)
        return rho * prop['R'] * T * (1 + prop['delta'] * prop['ar_delta'])
    
    elif param == 'u':  # Specific internal energy (J/kg)
        prop = explicit_helmholtz_derivs(rho, T)
        return prop['R'] * T * prop['tau'] * (prop['ao_tau'] + prop['ar_tau'])
    
    elif param == 's':  # Specific entropy (J/kg*K)
        prop = explicit_helmholtz_derivs(rho, T)
        return prop['R'] * (prop['tau'] * (prop['ao_tau'] + prop['ar_tau']) - prop['ao'] - prop['ar'])
    
    elif param == 'h':  # Specific enthalpy (J/kg)
        prop = explicit_helmholtz_derivs(rho, T)
        return prop['R'] * T * (1 + prop['tau'] * (prop['ao_tau'] + prop['ar_tau']) + prop['delta'] * prop['ar_delta'])
    
    elif param == 'a':  # Helmholtz free energy per unit mass (J/kg)
        prop = explicit_helmholtz_derivs(rho, T)
        return prop['R'] * T * (prop['ao'] + prop['ar'])
    
    elif param == 'mu': # Chemical Potential (J/kg)
        P = thermo_span_wagner(rho=rho, T=T, param='p')
        a = thermo_span_wagner(rho=rho, T=T, param='a')
        return a + P / rho
    
    elif param == 'cv':  # Specific heat constant volume (J/kg*K)
        prop = explicit_helmholtz_derivs(rho, T)
        return prop['R'] * -prop['tau']**2 * (prop['ao_tautau'] + prop['ar_tautau'])
    
    elif param == 'cp':  # Specific heat constant pressure (J/kg*K)
        prop = explicit_helmholtz_derivs(rho, T)
        return prop['R'] * (-prop['tau']**2 * (prop['ao_tautau'] + prop['ar_tautau']) + (1 + prop['delta'] * prop['ar_delta'] - prop['delta'] * prop['tau'] * prop['ar_deltatau'])**2 / (1 + 2 * prop['delta'] * prop['ar_delta'] + prop['delta']**2 * prop['ar_deltadelta']))
    
    elif param == 'du_drho_const_T': # 
        prop = explicit_helmholtz_derivs(rho, T)
        return prop['R'] * T / rho * ( prop['tau'] * prop['delta'] * prop['ar_deltatau'])
    
    elif param == 'dP_dT_const_rho':
        prop = explicit_helmholtz_derivs(rho, T)
        return rho * prop['R'] *(1 + prop['delta'] * prop['ar_delta'] - (prop['tau'] * prop['delta'] * prop['ar_deltatau']) )
    
    elif param == 'dP_drho_const_T':
        prop = explicit_helmholtz_derivs(rho, T)
        return prop['R'] * T * (1 + 2 * prop['delta'] * prop['ar_delta'] + (prop['delta']**2) * prop['ar_deltadelta'] )
    
    elif param == 'd_rho_dT_const_P':
        prop = explicit_helmholtz_derivs(rho, T)
        return (rho * prop['R'] *(1 + prop['delta'] * prop['ar_delta'] - (prop['tau'] * prop['delta'] *prop['ar_deltatau']) )) / (prop['R'] * T * (1 + 2 * prop['delta'] * prop['ar_delta'] + (prop['delta']**2) * prop['ar_deltadelta'] ))
    
    

    elif param == 'P_sat':
        pcnt = 0.1
        P_sat_est = P_sat_anc(T) # Start: Ancillary Eqn to sol an estimate for saturation temp
        while pcnt < 0.5:
            try:
                P_min, P_max = sol_init_guess(P_sat_est, pcnt)
                return brentq(lambda P: span2000_residual(P, T), P_min, P_max, xtol=1e-6, rtol=1e-6, maxiter=100)
            except ValueError:
                pcnt +=0.05
        raise RuntimeError("Root finding failed")
    
    elif param == "T_sat":
        T_guess = T_sat_anc(P_input)        
        return secant(lambda T: thermo_span_wagner(T=T, param="P_sat") - P_input, T_guess)

    elif param == "h_sat_vap":
        roots = find_all_sweos_density_roots(T, P, rho_min=0.1, rho_max=1200,n_points=1000)
        return thermo_span_wagner(rho=roots[0], T=T, param='h') # roots[0] Lowest density root (vapor)

    elif param == "h_sat_liq":
        roots = find_all_sweos_density_roots(T, P, rho_min=0.1, rho_max=1200,n_points=1000)
        return thermo_span_wagner(rho=roots[-1], T=T, param='h') # roots[-1] Highest density root (liquid)

    else:
        raise NotImplementedError(f'{param} is not implemented or incorrectly entered, see thermo_span_wagner()')



















"""T_input = 293
sol_psat = thermo_span_wagner(T=T_input, param="P_sat")
print("tsw P_sat:", sol_psat)"""

P_input = 4e6


sol_tsat = thermo_span_wagner(P=P_input, param="T_sat")
sol_hsat_liq = thermo_span_wagner(P=P_input, T=sol_tsat, param="h_sat_liq")
sol_hsat_vap = thermo_span_wagner(P=P_input, T=sol_tsat, param="h_sat_vap")
print("tsw T_sat:", sol_tsat, sol_hsat_liq, sol_hsat_vap)




"""
import CoolProp.CoolProp as CP
print("coolprop P_sat: ", CP.PropsSI('P', 'T', T_input, 'Q', 0, 'N2O'))

import numpy as np
import matplotlib.pyplot as plt

# Inputs
P_sat = thermo_span_wagner(T=T_input, param="P_sat")  # [Pa], your given saturation pressure

# Density sweep range
rhos = np.linspace(0.1, 1200, 1000)  # [kg/m^3]

mu_vals = []
for rho in rhos:
    try:
        mu_vals.append( thermo_span_wagner(rho=rho, T=T_input, param="mu") )
    except Exception as e:
        mu_vals.append(np.nan)

mu_vals = np.array(mu_vals)

# Plot
plt.figure(figsize=(10,6))
plt.plot(rhos, mu_vals, label='mu vs rho')
plt.xlabel('Density [kg/m³]')
plt.ylabel('Chemical Potential [J/kg]')
plt.title(f'Chemical Potential at T = {T_input} K, Sweeping Density')
plt.grid(True)

# Highlight regions where P matches P_sat (roots of EOS)
P_residuals = [thermo_span_wagner(rho=rho, T=T_input, param='p') - P_sat for rho in rhos]
signs = np.sign(P_residuals)
crossings = np.where(np.diff(signs))[0]

for idx in crossings:
    plt.axvline(rhos[idx], color='red', linestyle='--', alpha=0.6)
    print("root: ", rhos[idx])

plt.legend()
plt.show()


rhos = np.linspace(50, 1200, 300)

shift = thermo_span_wagner(rho=rho, T=T_input, param="a") - CP.PropsSI('HELMHOLTZMASS', 'D', rho, 'T', T_input, 'N2O')

a_model = np.array([thermo_span_wagner(rho=rho, T=T_input, param="a") for rho in rhos])
a_cp = np.array([CP.PropsSI('HELMHOLTZMASS', 'D', rho, 'T', T_input, 'N2O') + shift for rho in rhos])

# Plot
plt.plot(rhos, a_cp, label='CoolProp', linewidth=2)
plt.plot(rhos, a_model, '--', label='tswa', linewidth=2)
plt.xlabel("Density [kg/m³]")
plt.ylabel("Helmholtz Energy a [J/kg]")
plt.title(f"Helmholtz Energy vs Density at T = {T_input} K")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
"""