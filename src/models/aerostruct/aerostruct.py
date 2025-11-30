import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
from uvicrocketpy.mathutils.vector_matrix import Vector, Matrix


# NOTE: DOES NOT WORK! BAD

@dataclass
class PointLoad:
    force: float
    r: float



def build_inertial_loads(t, mass_data, rocket, flight):

    ## setup

    CG = rocket.center_of_mass(t)

    print("CG", CG)

    a_normal = flight.aerodynamic_lift(t)/rocket.total_mass(t)
    

    # sol point mass inertia along 1D rkt #TODO: INCLUDE PROPELLANT MASS!!!
    Ic_total = 0
    for m in mass_data:
        r_k = m.cg[-1] - CG
        Ic_total += m.mass*r_k**2

    angular_accel = flight.aerodynamic_spin_moment(t) /Ic_total


    print("a_normal: ", flight.aerodynamic_lift(t),rocket.total_mass(t) )
    print("angular_accel: ", angular_accel, flight.aerodynamic_spin_moment(t), Ic_total)

    inertial_loads = []

    for  m in mass_data:
        r_k = m.cg[-1] - CG
        a_rotate = -r_k*angular_accel
        a_k = -a_normal+a_rotate

        F = m.mass * a_k
        inertial_loads.append( PointLoad(force=F, r=m.cg[-1]) )

        print("look at force: ", m.mass, a_k, a_rotate, a_normal, F )        

    return inertial_loads

def build_aero_loads(t, rocket, flight):

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
        F = np.sqrt(X**2 + Y**2)

        aero_loads.append( PointLoad(force=F, r=(CG-comp_cp[-1])) )

    return aero_loads

def lump_point_loads(point_loads, tol=1e-6):
    """
    Combine point loads with identical (or nearly identical) r positions.
    Returns a new list of PointLoad objects.
    """
    # Sort by position
    point_loads = sorted(point_loads, key=lambda p: p.r)

    lumped = []
    current_r = None
    current_force = 0.0

    for p in point_loads:
        if current_r is None:
            current_r = p.r
            current_force = p.force
            continue

        # If positions are "the same" → lump them
        if abs(p.r - current_r) < tol:
            current_force += p.force
        else:
            lumped.append(PointLoad(force=current_force, r=current_r))
            current_r = p.r
            current_force = p.force

    # add last group
    lumped.append(PointLoad(force=current_force, r=current_r))

    return lumped


def sol_shear_n_bending(CG, point_loads, beam_length):

    # sol boundary conditions
    point_loads.append(PointLoad(force=0.0, r=0.0))

    total_force = 0.0
    for p in point_loads:
        total_force+=p.force

    point_loads.append(PointLoad(force=-total_force, r=CG))


    # sort by arm position
    point_loads = lump_point_loads(point_loads)

    xs =[]
    V = []
    M = []

    V_x = 0.0
    r_prev = 0.0
    M_x = 0.0
    for p in point_loads:
        xs.append(p.r)

        V_x -= p.force
        V.append(V_x)

        M_x +=V_x * (p.r-r_prev)
        M.append(M_x)

        r_prev = p.r

    return xs, V, M

def plot_shear_n_bending(xs, V, M):
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





def aerostruct(rkt_length, mass_data, rocket, flight):


    t_maxq = flight.max_dynamic_pressure_time

    ### build inertial loads
    dry_inertial_loads = build_inertial_loads(t_maxq, mass_data, rocket, flight)


    ### build aero loads
    aero_loads = build_aero_loads(t_maxq, rocket, flight)

    ### plot/summarize shear n bending from point loads
    point_loads = dry_inertial_loads + aero_loads
    CG = rocket.center_of_mass(t_maxq)

    print("aoa: ", flight.angle_of_attack(t_maxq))

    x, V, M = sol_shear_n_bending(CG, point_loads, rkt_length)
    plot_shear_n_bending(x, V, M)



    return 0