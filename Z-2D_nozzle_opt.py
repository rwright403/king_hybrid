import numpy as np
import matplotlib.pyplot as plt


### inputs: 
elevation = 285 #m
apogee_height = 6096 #m
optimal_height = 1 #m

# MAX_T_AMB
P_cc = 28e5      # chamber pressure [Pa]

# MIN_T_AMB
#P_cc = 23.5e5

y = 1.23                # ratio of specific heats

P_exit = 1000 * 101.29*( ((15.04 - 0.00649*optimal_height)+273.1)/288.08)**5.256 #Pa
Cf_opt = np.sqrt( ((2*y**2)/(y-1)) * ( (2/(y+1))**((y+1)/(y-1)) ) * (1- (P_exit/P_cc)**((y-1)/y)) ) 

# === Step 3: Evaluate Nozzle Performance  (OPTIMIZED DESIGN ALTITUDE) ===
# keep your existing A_throat sizing from It_est / tburn / Cf_opt
A_t = 0.25*np.pi*( 0.036673386 )**2      #It_est / (tburn * P_cc * Cf_opt)
throat_diam = 39.3701 * (2.0 * np.sqrt(A_t/np.pi))


# Altitude grid for the area calculation
h0 = float(elevation)
H  = float(elevation + apogee_height)
height_arr = np.linspace(h0, H, 300)

# --- local helpers (fine to be nested) ---
def P_atm(h):
    """ISA-ish model (Pa)."""
    t = 15.04 - 0.00649*h  # Celsius
    return 1000.0 * 101.29 * ((t + 273.1) / 288.08)**5.256

def Cf_isent(P_exit_local):
    """Isentropic Cf for given exit static pressure (Pa)."""
    pr = np.clip(P_exit_local / P_cc, 1e-12, 0.999999)
    term = 1.0 - pr**((y - 1.0)/y)
    return np.sqrt((2*y**2/(y-1.0)) * (2/(y+1.0))**((y+1.0)/(y-1.0)) * term)

def eps_from_Pexit(P_exit_local):
    """Expansion ratio eps = A_e/A_t for given P_exit and P_cc."""
    pr = np.clip(P_exit_local / P_cc, 1e-12, 0.999999)
    num = ((y+1.0)/2.0)**(1.0/(y-1.0)) * pr**(1.0/y)
    den = np.sqrt(((y+1.0)/(y-1.0)) * (1.0 - pr**((y-1.0)/y)))
    return 1.0 / (num * den)

def build_curves(h_star):
    """
    Return thrust arrays for:
    - blue: fixed nozzle sized at design altitude h_star (P_exit_star = P_atm(h_star))
    - orange: always-optimally-expanded (P_exit = P_atm(h))
    Also returns eps(h_star).
    """
    P_exit_star = P_atm(h_star)         # design exit pressure (fixed nozzle)
    eps_star    = eps_from_Pexit(P_exit_star)
    A_e         = eps_star * A_t

    P_a = P_atm(height_arr)

    # Orange: always optimally expanded (P_exit = P_a)
    Cf_opt_all = Cf_isent(P_a)

    # Blue: single nozzle sized at h_star (add pressure mismatch term)
    Cf_fixed = Cf_isent(P_exit_star) + (P_exit_star - P_a)/P_cc * eps_star

    T_opt   = P_cc * A_t * Cf_opt_all
    T_fixed = P_cc * A_t * Cf_fixed
    return T_fixed, T_opt, eps_star, P_exit_star

def objective(h_star, L2=False):
    """Area between curves over altitude (default L1)."""
    T_fixed, T_opt, _, _ = build_curves(h_star)
    diff = T_fixed - T_opt
    if L2:
        return np.trapz(diff**2, height_arr)
    else:
        return np.trapz(np.abs(diff), height_arr)

# --- optimize design altitude h* to minimize area ---
try:
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(lambda h: objective(h, L2=False),
                        bounds=(h0, H), method="bounded",
                        options={"xatol": 1.0})
    h_star_opt = float(res.x)
    area_min   = float(res.fun)
except Exception:
    # Fallback: grid search if SciPy isn't available
    candidates = np.linspace(h0, H, 401)
    vals = np.array([objective(h) for h in candidates])
    idx = int(np.argmin(vals))
    h_star_opt = float(candidates[idx])
    area_min   = float(vals[idx])

# Build optimum curves and parameters
T_fixed_opt, T_opt_opt, eps_opt, P_exit_star = build_curves(h_star_opt)
A_exit = eps_opt * A_t
expratio = eps_opt  # use optimized expansion ratio downstream

#update throat diam:
throat_diam = np.sqrt(A_t/np.pi)

# Plot vertical line for minimum T/W
#min_start_thrust = (rocket_dry_mass*9.81) * min_TW_ratio

print("\n=== Nozzle design-altitude optimization ===")
print(f"Optimal design altitude h*: {h_star_opt:.1f} m")
print(f"Design exit pressure P_exit(h*): {P_exit_star:.0f} Pa")
print(f"Expansion ratio eps = A_e/A_t: {eps_opt:.3f}")
print(f"Throat Area {throat_diam:.6f} m, [{39.3701*throat_diam:.6f}] in corresponding to Exit area A_exit: {A_exit:.6f} m^2")
print(f"Area between curves (L1): {area_min:.3e} N·m")
print(f"Average mismatch per meter: {area_min/(H-h0):.2f} N")

# Plot (blue = fixed @ h*, orange = always-optimal)
plt.figure()
plt.plot(T_fixed_opt, height_arr, label=f'fixed nozzle @ h*={h_star_opt:.0f} m')
plt.plot(T_opt_opt,   height_arr, label='optimal expansion throughout burn')
#plt.axvline(x=min_start_thrust, color='r', label=f'min starting thrust for T/W of {min_TW_ratio}')
plt.title(f'Altitude vs Thrust (throat_diam {throat_diam:.3f}\" ; eps={eps_opt:.2f})')
plt.xlabel('Thrust (N)'); plt.ylabel('Altitude (m)')
plt.grid(True); plt.legend(); plt.show()