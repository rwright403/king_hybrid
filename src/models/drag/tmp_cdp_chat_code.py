
import numpy as np

def _tangent_ogive_profile(L, R, x):
    """
    Tangent ogive radius r(x) and slope dr/dx for 0<=x<=L (to cylinder radius R).
    """
    Rs = (R**2 + L**2) / (2.0*R)       # sphere radius
    xc = L - Rs                        # sphere center x
    under = np.clip(Rs**2 - (x - xc)**2, 1e-30, None)
    r = np.sqrt(under) - (Rs - R)
    drdx = -(x - xc)/np.sqrt(under)
    return r, drdx
 
def _cone_pressure_linearized(M, gamma, p_inf, delta):
    """
    Very small-angle *axisymmetric* cone surface pressure.
    Linearized supersonic estimate:  Cp ≈ (2/β) * sin(delta)
    p = p_inf + Cp * q_inf,  q_inf = 0.5*γ*M^2*p_inf
    (Feel free to swap with your Taylor–Maccoll cone solver.)
    """

    beta = np.sqrt(M**2 - 1.0)
    Cp = (2.0/beta) * np.sin(delta)
    q_inf = 0.5*gamma*M**2*p_inf
    return p_inf + Cp*q_inf

    
def cdp_supersonic_shockexp(M, gamma, L_nose, R, npts=200, p_inf=101325.0, cone_pressure_fn=_cone_pressure_linearized):
    """
    Supersonic (M>1) **pressure/wave drag** for an axisymmetric body nose
    using the 2nd-order shock-expansion *marching* method.

    Geometry: tangent-ogive nose of length L_nose to cylinder radius R.
    Reference area A_ref = πR^2 (Barrowman convention).

    Parameters
    ----------
    M, gamma : float
        Freestream Mach and specific-heat ratio.
    L_nose, R : float
        Nose length and reference radius (D/2).
    npts : int
        Number of axial stations along the nose (panels = npts-1).
    p_inf : float
        Freestream static pressure [Pa].
    cone_pressure_fn : callable(M, gamma, p_inf, delta)->p_c
        Returns cone surface pressure for a cone with slope `delta`.

    Returns
    -------
    C_Dp : float
        Zero-lift pressure (wave) drag coefficient referenced to A_ref=πR^2.
    extras : dict
        Useful arrays for debugging/plots.
    """

    # Discretize nose
    x = np.linspace(0.0, L_nose, npts)
    # avoid singular tip in slope evaluation
    x[0] = 1e-9*L_nose
    r, drdx = _tangent_ogive_profile(L_nose, R, x)
    delta = np.arctan(drdx)             # surface angle to axis
    dx = np.diff(x)

    # Constants in the marching formulas
    beta = np.sqrt(M**2 - 1.0)
    B = gamma*p_inf*M**2/(2.0*(M**2 - 1.0))       # same at all stations for uniform freestream
    Omega = (1.0/M) * (1.0 + 0.5*(gamma-1.0)*M**2)**((gamma+1.0)/(2.0*(gamma-1.0)))

    # Allocate arrays (edge indexing: 1=prev rear, 2=current front, 3=current rear)
    p2 = np.empty(npts-1)       # surface pressure at front edge of each panel
    p3 = np.empty(npts-1)       # surface pressure at rear  edge of each panel
    dpds_1 = 0.0                # incoming from "previous" (initialize at tip)
    d1 = delta[0]               # previous trailing slope; at tip, use first slope

    # Initialize front-edge pressure at the tip: start from freestream
    # (If you have a better tip condition, replace here.)
    p_front = p_inf

    # Accumulator for integral
    CDP_sum = 0.0

    for j in range(npts-1):
        # panel ends
        x2, x3 = x[j], x[j+1]
        r2, r3 = r[j], r[j+1]
        d2, d3 = delta[j], delta[j+1]
        dxj = dx[j]

        # "cone" pressure evaluated with the local front slope
        p_c2 = cone_pressure_fn(M, gamma, p_inf, d2)

        # ---- (A) front-edge gradient (Eq. 3-70) ----
        # with uniform freestream B1=B2=B, Ω1=Ω2=Ω; keep symbols for clarity
        dpds_2 = (B/r2)*((Omega/Omega)*np.sin(d1) - np.sin(d2)) + (B/B)*(Omega/Omega)*dpds_1

        # ---- pressure field along the panel (Eq. 3-68) ----
        # exponential coordinate η
        denom = (p_c2 - p_front) * np.cos(d2)
        # guard against divide-by-zero if someone supplies p_front≈p_c2
        denom = denom if abs(denom) > 1e-24 else np.sign(denom)*1e-24

        # midpoint for integration
        xm = 0.5*(x2 + x3)
        eta_mid = dpds_2 * (xm - x2) / denom
        p_mid = p_inf - (p_c2 - p_front) * np.exp(-eta_mid)

        # ---- accumulate drag integral (axisymmetric, φ-integral→π) ----
        P_hat_mid = (p_mid - p_inf)/p_inf
        CDP_sum += P_hat_mid * np.sin(d2) * dxj

        # ---- (rear-edge values to carry) ----
        eta_rear = dpds_2 * (x3 - x2) / denom
        p_rear = p_inf - (p_c2 - p_front) * np.exp(-eta_rear)

        # gradient carry-over (Eq. 3-71)
        dpds_3 = ((p_c2 - p_rear)/(p_c2 - p_front)) * dpds_2

        # record and advance
        p2[j] = p_front
        p3[j] = p_rear
        dpds_1 = dpds_3   # becomes incoming for next panel
        d1 = d3           # current trailing slope is next "previous" slope
        p_front = p_rear  # next panel front-edge pressure

    # Final coefficient
    C_Dp = (4.0*np.pi/(gamma*M**2)) * CDP_sum

    #extras = dict(x=x, r=r, delta=delta, p2=p2, p3=p3)
    return C_Dp #, extras
