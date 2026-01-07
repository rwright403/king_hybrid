import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
from uvicrocketpy.mathutils.vector_matrix import Vector, Matrix



@dataclass
class PointLoad:
    force: float
    r: float

def rocketpy_load_case(t, mass_data, rocket, flight):
    #DO NOT USE, DOES NOT WORK WELL angular accel too small

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

    print("check len mass data and inertial loads:", len(mass_data), len(inertial_loads) )     

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

def build_aero_loads_openrocket(flight):

    #LIST OF CN_a from OPENROCKET:
    or_aero_comp = [ ORCompAnalysis(x_cp=4.26305, cna=2.286 ), # nosecone
                    ORCompAnalysis(x_cp=1.98 , cna=3.94 ), # body tube
                    ORCompAnalysis(x_cp=0.2568 , cna=6.354), # fins
    
    ]
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

        print(f"shear n bend: {(p.r):.3f}, {V_x:.3f}, {M_x:.3f},    {p.force:.3f}")

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





def aerostruct(rkt_length, mass_data, rocket, flight):

    t_maxq = flight.max_dynamic_pressure_time
    CG = rocket.center_of_mass(t_maxq) #does this include propellants?

    
    angular_accel, a_normal = rocketpy_load_case(t_maxq, mass_data, rocket, flight)



    ### build aero loads
    #aero_loads = build_aero_loads_rocketpy(t_maxq, rocket, flight)
    aero_loads = build_aero_loads_openrocket(flight)

    aoa =flight.angle_of_attack(t_maxq)
    print("aoa: ", aoa)


    CG = rocket.center_of_mass(t_maxq) #does this include propellants?
    print("CG: ", CG, rkt_length)

    angular_accel, a_normal = or_load_case(t_maxq, CG, mass_data, rocket, aero_loads)


    ### build inertial loads
    dry_inertial_loads = build_inertial_loads(CG, mass_data, angular_accel, a_normal)

    x, V, M = sol_shear_n_bending(CG, dry_inertial_loads, aero_loads, rkt_length)
    plot_shear_n_bending(x, V, M)
    
