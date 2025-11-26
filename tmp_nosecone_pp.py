import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.interpolate import RegularGridInterpolator

# ----------------------------
# Geometry helper
# ----------------------------
def _tangent_ogive_profile(L, R, x):
    """
    Return r(x) and dr/dx for a tangent ogive of length L and radius R.
    """
    Rc = (L**2 + R**2) / (2*R)     # ogive sphere radius
    under = Rc - R

    h = np.sqrt(Rc**2 - (x - L)**2)
    r = h - under

    drdx = np.zeros_like(r)
    drdx[:-1] = (r[1:] - r[:-1]) / (x[1:] - x[:-1])
    drdx[-1] = 0.0
    return r, drdx




angles = np.array([2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0])         # degrees
M_inf_vals = np.array([1.5, 1.75, 2.0, 2.5])
# ----------------------------
# Table 3-1 Initial Cone Surface
# ----------------------------

# localMach[i][j] = M_local at M_inf_vals[i], angles[j]
localMachCones = np.array([
#  delta = 2.5,         5.0,      7.5,     10.0,        12.5,     15.0,      17.5,      20.0,      22.5,       25.0,    27.5,       30.0
    [   1.4866884, 1.4579318, 1.4197316, 1.3748579, 1.3248992, 1.2707392, 1.2127162, 1.1506488, 1.0837313, 1.0100703, 1.0100703, 1.0100703],   # 1.5
    [   1.7340335, 1.7004643, 1.6568043, 1.6064462, 1.5513276, 1.4925484, 1.4306418, 1.3657361, 1.2976505, 1.2259024, 1.1495789, 1.0669388],   # 1.75
    [   1.9808615, 1.9415254, 1.8912335, 1.8340453, 1.7721933, 1.7068837, 1.6386764, 1.5677572, 1.4940884, 1.4174680, 1.3375202, 1.2535953],   # 2.0
    [   2.4729472, 2.4193976, 2.3528628, 2.2788768, 2.2001554, 2.1179422, 2.0327606, 1.9448503, 1.8543543, 1.7613748, 1.6659616, 1.5680718],   # 2.5
])

# Build interpolator ONCE (don't rebuild inside lookup function)
interpLocalMachCones = RegularGridInterpolator(
    (M_inf_vals, angles),
    localMachCones,
    bounds_error=False,
    fill_value=None
)

# localMach[i][j] = M_local at M_inf_vals[i], angles[j]
surfPressureCoeffCones = np.array([
#  delta = 2.5,         5.0,      7.5,     10.0,        12.5,     15.0,      17.5,      20.0,      22.5,       25.0,    27.5,       30.0
    [   0.0123382, 0.0396615, 0.0773864, 0.1237983, 0.1781167, 0.2400313, 0.3095450, 0.3869791, 0.4731538, 0.5699888, 0.5699888, 0.5699888],   # 1.5
    [   0.0114414, 0.0363287, 0.0704289, 0.1122701, 0.1612014, 0.2169225, 0.2792900, 0.3482459, 0.4238192, 0.50619022, 0.59585744, 0.69405243],   # 1.75
    [   0.0107862, 0.0339457, 0.0655812, 0.1044583, 0.1500854, 0.2022339, 0.2607607, 0.3255312, 0.39641336, 0.47330179, 0.55616008, 0.64511954],   # 2.0
    [   0.00982981, 0.0305876, 0.05895952, 0.09413081, 0.1358608, 0.1840343, 0.23852578, 0.29915906, 0.36571047, 0.43791904, 0.51550425, 0.59818912],   # 2.5
])

# Build interpolator ONCE (don't rebuild inside lookup function)
interpSurfPressureCoeffCones = RegularGridInterpolator(
    (M_inf_vals, angles),
    surfPressureCoeffCones,
    bounds_error=False,
    fill_value=None
)


# ----------------------------
# table_lookup function
# ----------------------------
def table_lookup(Mach, delta_rad, interp):
    delta_deg = np.degrees(delta_rad)

    # Interpolator requires a 2D input (N,2)
    val = interp([[Mach,delta_deg]])[0]

    return float(val)




# ----------------------------
# beta and omega functions (from eq. 3-70 notation)
# ----------------------------
def beta_coeff(M, gamma):
    """
    beta from Barrowman eq. (3-70):
    """
    return (gamma*M**2) / (2*(M**2 - 1.0))


def omega_coeff(M, gamma):
    """
    omega from Barrowman eq. (3-70):
    """
    return (1.0/M) * ((1.0 + 0.5*(gamma-1.0)*M**2) / (0.5*(gamma+1.0)) )**((gamma+1.0)/(2.0*(gamma-1.0)))




# ================================================================
# FINAL FUNCTION: 2nd-order Shock-Expansion marching CDp
# ================================================================
def shockexp_pressure_distribution(M_inf, p_inf, L_nose, R, gamma=1.4, npts=80):

    # Discretize nose
    theta = np.linspace(0, np.pi/2, npts)
    x = L_nose * (1 - np.cos(theta))


    r, drdx = _tangent_ogive_profile(L_nose, R, x)
    delta = np.arctan(drdx)

    # Output arrays
    dpds_front = np.zeros(npts-1)
    p2_arr = np.zeros(npts-1)
    p3_arr = np.zeros(npts-1)
    pmid_arr = np.zeros(npts-1)

    # freestream dynamic pressure:
    q_inf = 0.5*gamma*(M_inf**2)*p_inf

    j = 0 

    d1 = delta[j]                            # slope of panel 1
    Cp1 = table_lookup(M_inf, d1, interpSurfPressureCoeffCones)
    p_front = p_inf + Cp1 * q_inf              # p2 for panel 1

    # Barrowman boundary condition:
    dpds_1 = 0.0

    j+=1 # skip first panel

    # March from j = 1 to end
    for j in range(j, npts-1):

        x2, x3 = x[j], x[j+1]
        r2 = r[j]
        d2, d3 = delta[j], delta[j]


        # Local Mach numbers
        M1 = table_lookup(M_inf, d1, interpLocalMachCones)
        M2 = table_lookup(M_inf, d2, interpLocalMachCones)

        # Beta and Omega
        B1 = beta_coeff(M1, gamma)
        B2 = beta_coeff(M2, gamma)
        omega1 = omega_coeff(M1, gamma)
        omega2 = omega_coeff(M2, gamma)

        # Cone pressure at front
        Cp = table_lookup(M_inf, d2, interpSurfPressureCoeffCones)
        p_c2 = p_inf + Cp * q_inf

        #print("p_c2: ", p_c2, p_inf, Cp, q_inf, d2)

        # Eq (3-70) - dp/ds at front
        dpds_2 = (B2/r2)*(omega1*np.sin(d1) - omega2*np.sin(d2)) + (B2/B1)*(omega1/omega2)*dpds_1

        #print(omega1, omega2, B1, B2, M1, M2, M_inf)

        #print(d1, d2)

        #print( (B2/r2), omega1*np.sin(d1), -omega2*np.sin(d2)  )

        #print(f"{(B2/r2):.3f}, {(omega1*np.sin(d1) - omega2*np.sin(d2)):.3f}, {(B2/B1)*(omega1/omega2)*dpds_1}")
        dpds_front[j] = dpds_2

        # Eq (3-68) - pressure along panel
        xm = 0.5*(x2 + x3)
        eta_mid = dpds_2 * (xm - x2)/ ((p_c2 - p_front)*np.cos(d2))

        #print( (p_c2 - p_front)*np.cos(d2), p_c2,p_front )



        p_mid = p_front + (p_c2 - p_front)*np.exp(-eta_mid)
        pmid_arr[j] = p_mid

        eta_rear = dpds_2 * (x3 - x2)/ ((p_c2 - p_front)*np.cos(d2))
        #print(f"{eta_rear:.10f}, {dpds_2:.3f}, {(x3-x2):.3f}, {((p_c2 - p_front)*np.cos(d2))}, {np.cos(d2)}")

        #print(p_c2, p_front)



        p_rear = p_front + (p_c2 - p_front)*np.exp(-eta_rear)

        print(p_rear, p_front, - (p_c2 - p_front)*np.exp(-eta_rear), (p_c2 - p_front), eta_rear)




        p2_arr[j] = p_front
        p3_arr[j] = p_rear

        # Eq (3-71) - slope carryover
        dpds_3 = ((p_c2 - p_rear)/(p_c2 - p_front))*dpds_2

        #print(f"{j:.3f}, {d1:.3f}, {d2:.3f}, {np.sin(d1):.3f}, {np.sin(d2):.3f}, {omega1:.3f}, {omega2:.3f}, {dpds_2:.3f}")

        #print(f"{j:.3f}, {p_mid:.0f}, {p_inf:.0f}, {-(p_c2 - p_front):.0f}, {p_c2:.0f}, {p_front:.0f}, {np.exp(-eta_rear):.0f}")

        # Advance marching state
        dpds_1 = dpds_3
        p_front = p_rear
        d1 = d3


    return {
        "x": x,
        "r": r,
        "delta": delta,
        "dpds_front": dpds_front,
        "p2": p2_arr,
        "p3": p3_arr,
        "pmid": pmid_arr
    }





def plot_shockexp_pressure_contour(data, field="pmid"):
    """
    Plot the tangent-ogive contour r(x) colored by surface pressure.

    field options:
        "pmid"  -> midpoint panel pressure (recommended)
        "p2"    -> front-edge pressure
        "p3"    -> rear-edge pressure
    """

    x = data["x"]
    r = data["r"]      # radius distribution
    p = data["pmid"]

    # pressure is panel-based (N-1), need vertex-based (N)
    p_vertices = np.zeros_like(x)
    p_vertices[:-1] = p
    p_vertices[-1] = p[-1]

    # Build line segments for contour
    points = np.array([x, r]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Create colored line
    lc = LineCollection(segments, cmap="plasma")
    lc.set_array(p_vertices)
    lc.set_linewidth(3.0)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.add_collection(lc)

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(r.min(), r.max())
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel("x (m)")
    ax.set_ylabel("radius r(x) (m)")
    ax.set_title(f"Nosecone Surface Pressure Distribution")

    # Add colorbar
    plt.colorbar(lc, ax=ax, label="Surface Pressure (Pa)")

    plt.grid(True)
    plt.tight_layout()
    plt.show()





M_inf = 1.8
p_inf = 101325 
L_nose = (0.0254*20)
R = 0.5*(0.0254*5.5)

data = shockexp_pressure_distribution(M_inf, p_inf, L_nose, R)
plot_shockexp_pressure_contour(data)