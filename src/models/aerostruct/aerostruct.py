import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
from uvicrocketpy.mathutils.vector_matrix import Vector, Matrix



@dataclass
class PointLoad:
    force: float
    r: float

def rocketpy_load_case(t, mass_data, rocket, flight):

    CG = rocket.center_of_mass(t)

    print("CG", CG)

    a_normal = flight.aerodynamic_lift(t)/rocket.total_mass(t)
    

    # sol point mass inertia along 1D rkt #TODO: NOT INCLUDE PROPELLANT MASS!!!
    Ic_total = 0
    for m in mass_data:
        r_k = m.cg[0] - CG
        Ic_total += m.mass*r_k**2

    angular_accel = flight.aerodynamic_spin_moment(t) /Ic_total # (N m) / (kg m^2)
    
    print("a_normal: ", flight.aerodynamic_lift(t), rocket.total_mass(t) )
    print("angular_accel: ", angular_accel, flight.aerodynamic_spin_moment(t), Ic_total)

    return angular_accel, a_normal

def or_load_case(t, CG, mass_data, rocket, aero_loads):

    aero_lift = 0
    aero_moment = 0
    for p in aero_loads:
        aero_lift+=p.force
        aero_moment+=p.force*(p.r-CG)

    a_normal = aero_lift/rocket.total_mass(t)
    

    # sol point mass inertia along 1D rkt #TODO: INCLUDE PROPELLANT MASS!!!
    Ic_total = 0
    for m in mass_data:
        r_k = m.cg[0] - CG
        Ic_total += m.mass*r_k**2

    angular_accel =  aero_moment/Ic_total # (N m) / (kg m^2)
    

    return angular_accel, a_normal



def build_inertial_loads(CG, mass_data, angular_accel, a_normal):

    inertial_loads = []

    for m in mass_data:
        r_k = m.cg[0] - CG
        a_rotate = -r_k*angular_accel
        a_k = -a_normal+a_rotate

        F = m.mass * a_k
        inertial_loads.append( PointLoad(force=F, r=(r_k+CG) ) ) #NOTE: point loads in tail csys

        #print("look at force: ", m.mass, a_k, a_rotate, a_normal, F, r_k )     

    return inertial_loads



def build_aero_loads_rocketpy(t, rocket, flight):
    #returns aero loads about cg

    CG = rocket.center_of_mass(t)

    aero_loads = []
    # build rotation matrices with euler parameters:
    e = [flight.e0(t), flight.e1(t), flight.e2(t), flight.e3(t)]
    K = Matrix.transformation(e)
    Kt = K.transpose

    #NOTE: w in flight is: w = Vector([omega1, omega2, omega3])
        # check omega1 vs w1 - i think omega1 is the local var and w1 is the datastructure

    w = Vector([flight.w1(t),flight.w2(t),flight.w3(t)])
        
    # velocity_in_body_frame = Vector([vx_b, vy_b, vz_b])
    velocity_inertial_frame = Vector([flight.vx(t), flight.vy(t), flight.vz(t)])
    velocity_in_body_frame = Kt @ velocity_inertial_frame


    for aero_surface, _ in rocket.aerodynamic_surfaces:


        # Component cp relative to CDM in body frame
        comp_cp = rocket.surfaces_cp_to_cdm[aero_surface]
        # Component absolute velocity in body frame
        comp_vb = velocity_in_body_frame + (w ^ comp_cp)


        print("aero comp: ", comp_cp, comp_vb)

        # Wind velocity at component altitude
        comp_z = flight.z(t) + (K @ comp_cp).z
        comp_wind_vx = flight.env.wind_velocity_x.get_value_opt(comp_z)
        comp_wind_vy = flight.env.wind_velocity_y.get_value_opt(comp_z)
        # Component freestream velocity in body frame
        comp_wind_vb = Kt @ Vector([comp_wind_vx, comp_wind_vy, 0])
        comp_stream_velocity = comp_wind_vb - comp_vb
        comp_stream_speed = abs(comp_stream_velocity)
        comp_stream_mach = comp_stream_speed / flight.speed_of_sound(t)
        # Reynolds at component altitude
        # TODO: Reynolds is only used in generic surfaces. This calculation
        # should be moved to the surface class for efficiency
        comp_reynolds = (
            flight.env.density.get_value_opt(flight.z(t))
            * comp_stream_speed
            * aero_surface.reference_length
            / flight.env.dynamic_viscosity.get_value_opt(flight.z(t))
        )

        # Forces and moments
        X, Y, _, _, _, _ = aero_surface.compute_forces_and_moments(
            comp_stream_velocity,
            comp_stream_speed,
            comp_stream_mach,
            flight.density(t),
            comp_cp,
            w,
            comp_reynolds,
        )

        # resolve 3d aero loads to 1d aero load by taking magnitude of X Y
        F = X #np.sqrt(X**2 + Y**2)

        aero_loads.append( PointLoad(force=F, r=(comp_cp[-1]+CG)) )
        #print("aero load location: ",(CG-comp_cp[-1]), comp_cp[-1]+CG )

    return aero_loads

@dataclass
class ORCompAnalysis:
    x_cp: float
    cna: float

def build_aero_loads_openrocket():


    #LIST OF CN_a from OPENROCKET:
    """
    rkt_len = 4.15 #NOTE: OR in nose - tail csys, need to convert to our tail-nose
    or_aero_comp = [ ORCompAnalysis(x_cp=rkt_len-0.23621, cna=2.115 ), # nosecone
                    ORCompAnalysis(x_cp=rkt_len-2.325, cna=1.262 ), # body tube
                    ORCompAnalysis(x_cp=rkt_len-3.94235 , cna=6.277), # fins
    ]
    print("OR Coeffs from 2 deg aoa and 1.56 Ma. Rocket length is 4.15 [m]")
    print("NOTE: OR throws warning for body coeff in supersonic regime")
    """


    rkt_len = 4.15 #NOTE: OR in nose - tail csys, need to convert to our tail-nose
    or_aero_comp = [ ORCompAnalysis(x_cp=rkt_len-0.24124, cna=2.264 ), # nosecone
                    ORCompAnalysis(x_cp=rkt_len-2.325, cna=2.898 ), # body tube
                    ORCompAnalysis(x_cp=rkt_len-3.94235 , cna=6.524 ), # fins
    ]
    print("OR Coeffs from 4.6 deg aoa and 1.56 Ma. Rocket length is 4.15 [m]")
    print("NOTE: OR throws warning for body coeff in supersonic regime")
    


    L_ref = 0.141 #m
    A_ref = 0.0157 #m^2

    # 5 degree aoa taken at program's "worst wind dir". Mach number from MAXQ, rho_maxq and v_maxq as well

    rho_maxq = 0.989 #kg/m^3
    v_maxq = 425.498 #m/s

    aero_loads = []

    for comp in or_aero_comp:
        F_N = 0.5 *comp.cna* rho_maxq* (v_maxq**2)*A_ref
        aero_loads.append( PointLoad(force=F_N, r=comp.x_cp ))

    return aero_loads




def sol_shear_n_bending(CG, dry_inertial_loads, aero_loads, beam_length):

    point_loads = dry_inertial_loads + aero_loads

    # sol boundary conditions -csys 0 at cg
    point_loads.append(PointLoad(force=0.0, r=0 ))
    point_loads.append(PointLoad(force=0.0, r=beam_length ))

    total_force = 0.0
    for p in point_loads:
        total_force+=p.force

    point_loads.append(PointLoad(force= -total_force, r=CG))


    net_moment = 0.0
    for p in point_loads:
        net_moment-=p.force*(p.r-CG)

    # sort by arm position
    point_loads = sorted(point_loads, key=lambda p: p.r)

    xs = []
    V = []
    M = []

    V_x = 0.0
    r_prev = 0.0
    M_x = 0.0
    for p in point_loads:
        xs.append(p.r)

        #update moment first so V_x is previous
        if (p.r != CG):
            M_x += V_x * (p.r-r_prev)
        else: #moment at center of gravity equals the net moment about the cg from external aero loads (boundary condition)
            M_x += net_moment
            
        V_x -= p.force
        V.append(V_x)
        M.append(M_x)

        r_prev = p.r

    net_force = 0
    for p in point_loads:
        net_force+=p.force
    print("force: along beam should be 0: ", net_force)

    
    return xs, V, M



def plot_shear_n_bending(xs, V, M):
    """ plotting """
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    # --- Shear ---
    axes[0].plot(xs, V, linewidth=2)
    axes[0].axhline(0, color='black', linewidth=1)
    axes[0].set_ylabel("Shear Force V(x)")

    # --- Bending ---
    axes[1].plot(xs, M, linewidth=2)
    axes[1].axhline(0, color='black', linewidth=1)
    axes[1].set_ylabel("Bending Moment M(x)")
    axes[1].set_xlabel("Rocket Length (m)")

    plt.tight_layout()
    plt.show()


def aoa_from_velo_triangle(flight, t, v_gust = (150/3.6) ): # default [50 km/hr] converted to [m/s]

    aoa = np.degrees( np.arctan( v_gust/flight.free_stream_speed(t) ) )

    return aoa



def plot_point_loads(loads, beam_length):
    xs = [p.r for p in loads]
    Fs = [abs(p.force) for p in loads]

    maxF = max(Fs)
    scale = 0.25 / maxF

    fig, ax = plt.subplots()

    # Beam line
    ax.hlines(0, 0, beam_length, colors='black', linewidth=3)

    used_positions = []
    arrow_extents = []

    for p in loads:
        x = p.r
        arrow_len = p.force * scale

        # --- smart stagger logic ---
        level = 0
        min_spacing = 0.2
        while any(abs(x - ux) < min_spacing and level == l for ux, l in used_positions):
            level += 1
        used_positions.append((x, level))
        offset = level * 0.07
        # ---------------------------

        # arrow
        ax.annotate(
            "",
            xy=(x, arrow_len),
            xytext=(x, 0),
            arrowprops=dict(width=2, headwidth=10, headlength=12)
        )

        # label
        label_y = arrow_len * 1.15 + offset
        ax.text(
            x,
            label_y,
            f"{p.force:.1f} N",
            ha="center",
            va="bottom" if arrow_len > 0 else "top"
        )

        arrow_extents.append(label_y)

    # ---------- BIASED Y-LIMITS ----------
    max_y = max(arrow_extents + [0])
    min_y = min(arrow_extents + [0])

    # keep a tight band around beam, bias space upward
    lower_pad = 0.12          # small space below beam
    upper_pad = (max_y - min_y) * 0.35  # more space near labels

    ax.set_ylim(min_y - lower_pad, max_y + upper_pad)
    # ------------------------------------

    ax.set_xlim(-0.05 * beam_length, beam_length * 1.05)
    ax.set_xlabel("Beam Position (m)")
    ax.set_ylabel("Force Magnitude (scaled)")   # ← fixed 👍
    ax.grid(True)

    plt.show()


def print_point_loads(loads):
    print("\n")
    for p in loads:
        print(f"(location (tail to nose), force): {(p.r):.3f} [m], {p.force:.3f} [N]")
    print("\n")

def lookup_thrust_at_t(t, filepath):

    """
    Linearly interpolate y at x_query from a 2-column data file.
    First row is assumed to be a header and is skipped.
    """
    # load data, skipping header
    data = np.loadtxt(filepath, skiprows=1, delimiter=',')
    
    xs = data[:, 0]
    ys = data[:, 1]

    # ensure sorted by x
    sort_idx = np.argsort(xs)
    xs = xs[sort_idx]
    ys = ys[sort_idx]

    # numpy interp does linear interpolation
    # NOTE: outside range → it clamps to boundary values
    thrust = np.interp(t, xs, ys)
    print(f"Thrust at t = {t:.3f} [s]: Thrust = {thrust:.3f} [N] (careful, hardcoded filepath to results!!: {filepath})")





# Jack Parachute Differential Function
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


from scipy.integrate import solve_ivp

def make_reco_solver(m_parachute, mn, mr, v_wind):
    def sol_reco_ivp(rhoCdA, V):

        k = 4942  # shock cord spring const
        g = 9.81  # gravity
        t = np.linspace(0, 0.75, 5000)  # Time vector

        # Dynamic Parameters
        xp = 0                     # Parachute Initial Position            [m]     {always zero}
        xn = 0                     # Upper Airframe Initial Position       [m]     {always zero}
        xr = 0                     # Lower Airframe Initial Position       [m]     {always zero}
        xp_dot = 0                 # Parachute Initial Velocity            [m/s]   {zero unless reefed chute}
        xn_dot = V                 # Initial Velocity of Upper Airframe    [m/s]
        xr_dot = V                 # Initial Velocity of Lower Airframe    [m/s]

        y0 = [xp, xn, xr, xp_dot, xn_dot, xr_dot]


        solution = solve_ivp(
            dy_dt,
            (t[0], t[-1]),
            y0,
            args=(m_parachute, mn, mr, rhoCdA, k, g, v_wind),
            t_eval=t
        )

        # Extract results
        xp_sol = solution.y[0]
        xr_sol = solution.y[2] # want this one

        lower_shock_tension = k * (xr_sol - xp_sol)

        # largest magnitude
        max_tension = np.max(np.abs(lower_shock_tension))

        return max_tension

    return sol_reco_ivp


from tqdm import tqdm

def reco_load_map(flight, m_dry, mr=3.45, m_parachute=0.07, v_wind=(150/3.6)):

    mn = m_dry -(mr-m_parachute)

    rho_arr = []
    vz_arr = []
    
    t0 = flight.apogee_time
    t1 = flight.t
    n_points = 100

    # Create log-spaced values between 0 and 1, then scale to t0 --> t1
    #t_arr = t0 + (t1 - t0) * np.logspace(0, 1, n_points, base=10) / 10
    t_arr = np.linspace(t0, t1, n_points)

    vz_min = -150
    vz_max = -50

    for t in t_arr:
        vz = flight.vz(t)
        if (vz_min <= vz <= vz_max):
            rho_arr.append( flight.density(t) )
            vz_arr.append(vz)

        

    # Example inputs
    CdA_array = np.linspace(3, 152.7, 50) #np.linspace(0.3744, 152.7, 50)  # Cd * A values

    # Pre-allocate results
    force_map = np.zeros((len(vz_arr), len(CdA_array)))

    solver = make_reco_solver(m_parachute, mn, mr, v_wind)

    # Loop over all combinations
    for i in tqdm(range(len(vz_arr)), desc="Reco Combinations"):
        for j, CdA in enumerate(CdA_array):
        
            rho = rho_arr[i]
            vz = vz_arr[i]

            rhoCdA = rho*CdA

            # example masses and wind, replace with your actual values
            max_force = solver(rhoCdA, vz)
            force_map[i, j] = max_force

    # Make a color map
    plt.figure(figsize=(8,6))
    X, Y = np.meshgrid(CdA_array, vz_arr)
    plt.pcolormesh(X, Y, force_map, shading='auto', cmap='viridis')
    plt.colorbar(label='Max Lower Shock Tension [N]')
    plt.xlabel('CdA [m²]')
    plt.ylabel('Deployment Velocity [m/s]')
    plt.title('Max Lower Shock Tension vs CdA and Velocity')
    plt.show()










def aerostruct(rkt_length, mass_data, rocket, flight):


    ### MAX Q CASE!!! ###
    t_maxq = flight.max_dynamic_pressure_time
    CG = rocket.center_of_mass(t_maxq) #does this include propellants?

    
    angular_accel, a_normal = rocketpy_load_case(t_maxq, mass_data, rocket, flight)

    ### build aero loads
    #aero_loads = build_aero_loads_rocketpy(t_maxq, rocket, flight)
    aero_loads = build_aero_loads_openrocket()

    aoa = flight.angle_of_attack(t_maxq)
    print(f"nom flight aoa: , {aoa} degrees (for ref)")

    gust_aoa = aoa_from_velo_triangle(flight, t_maxq) # default [50 km/hr] converted to [m/s]
    print(f"\nOR Local Coeff inputs (aoa, Ma, neglect roll rate) : , {gust_aoa} degrees, {flight.mach_number(t_maxq)}\n")


    CG = rocket.center_of_mass(t_maxq) #does this include propellants?

    angular_accel, a_normal = or_load_case(t_maxq, CG, mass_data, rocket, aero_loads)


    ### build inertial loads
    dry_inertial_loads = build_inertial_loads(CG, mass_data, angular_accel, a_normal)

    ### plot loads:
    plot_point_loads( (dry_inertial_loads+aero_loads), rkt_length)

    x, V, M = sol_shear_n_bending(CG, dry_inertial_loads, aero_loads, rkt_length)

    print(f"Drag (aero force aligned w rocket axis) at t={t_maxq:.3f} [s], F_drag = {flight.aerodynamic_drag(t_maxq):.3f} [N]")
    lookup_thrust_at_t(t_maxq, "src/results/Skvader_config_35_FINAL_T_MAX/thrust.csv")
    plot_shear_n_bending(x, V, M)



    ### RAIL BURNOUT FORCE
    flight.plots.rail_buttons_forces()
    

    ### motor burnout case. t_burnout = 6.52 [s] for T_max thrust curve (worst case loading)
    t_burnout = 6.52
    print("\n\n")
    print(f"Drag (aero force aligned w rocket axis) at t={t_burnout:.3f} [s], F_drag = {flight.aerodynamic_drag(t_burnout):.3f} [N]")


    ### ESTIMATING RECO LOADS
    m_dry=0
    for m in mass_data:
        m_dry += m.mass
    reco_load_map(flight, m_dry)


