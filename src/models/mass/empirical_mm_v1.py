import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt


def diameter_from_volume_LD(V, L_over_D):
    D = (4*V / (np.pi * L_over_D))**(1/3)
    return D


@dataclass
class mass:
    mass: float
    cg: np.ndarray 
    inertia: np.ndarray # 3x3 Inertia Tensor


def hollow_cylinder_mass(rho, L, d_inner, d_outer):
    return rho * np.pi * L * ( (d_outer/2)**2 - (d_inner/2)**2 )

def hollow_cylinder_inertia(m, L, d_inner, d_outer):
    # Inertia of a hollow cylinder about itself    
    Ixx = 0.5 * m * ((d_outer/2)**2 + (d_inner/2)**2)
    Iyy = Izz = (1/12) * m * (3*((d_outer/2)**2 + (d_inner/2)**2) + L**2)
    inertia = np.array([[Ixx,0,0],[0,Iyy,0],[0,0,Izz]])
    return inertia


def hollow_cylinder_lumped_masses(rho, L, d_inner, d_outer, cyl_center, N=7):
    """
    Create N point-mass objects representing a hollow cylinder.

    Inputs:
        rho         : material density (kg/m^3)
        L           : length of cylinder (m)
        d_inner     : inner diameter (m)
        d_outer     : outer diameter (m)
        cyl_center  : np.array([x0, y0, z0]) global position of the cylinder's centroid
        N           : number of point masses

    Returns:
        list[mass]  : each mass has:
                        mass (kg)
                        cg (global coordinates)
                        inertia (about its own CG)
    """

    Ri = d_inner / 2
    Ro = d_outer / 2

    # Total cylinder mass
    total_mass = rho * np.pi * L * (Ro**2 - Ri**2)
    m_i = total_mass / N

    # Local positions along cylinder axis (centered at 0)
    xs_local = np.linspace(-L/2, L/2, N)

    # Local cross-section inertia (about *each lump's* CG)
    Ixx_local = 0.5 * m_i * (Ro**2 + Ri**2)
    Iyy_local = 0.25 * m_i * (Ro**2 + Ri**2)
    Izz_local = Iyy_local
    inertia_local = np.diag([Ixx_local, Iyy_local, Izz_local])

    mass_list = []
    for x_local in xs_local:
        cg_global = cyl_center + np.array([x_local, 0.0, 0.0])
        mass_obj = mass(
            mass=m_i,
            cg=cg_global,
            inertia=inertia_local.copy()   # still about the lump’s own CG
        )
        mass_list.append(mass_obj)

    return mass_list





def thin_cone_shell_mass(rho, H, d_outer, t=3e-3):
    r = 0.5 * d_outer
    s = np.sqrt(r**2 + H**2)      # slant height
    A = np.pi * r * s             # lateral surface area
    V = A * t                     # thin shell volume
    return rho * V
#NOTE: TEMP FOR DEBUGGING


def thin_hollow_cone_inertia(m, H, d_outer):
    R = d_outer / 2

    Ixx = 0.5 * m * R**2
    Iyy = 0.25 * m * R**2 + (1/12) * m * H**2
    Izz = Iyy
    inertia = np.array([[Ixx,0,0],[0,Iyy,0],[0,0,Izz]])
    return inertia



def total_fin_inertia(fin_span, fin_root_chord, fin_tip_chord, fin_thickness, fin_density, od_fuse, n_fins):

    
    # ------------------------------
    # 1. Single fin mass
    # ------------------------------
    c_avg = 0.5 * (fin_root_chord + fin_tip_chord)
    A_planform = fin_span * 0.5 * (fin_root_chord + fin_tip_chord)
    m = A_planform * fin_thickness * fin_density

    # local fin inertias (thin plate model)
    Ixx = (1/12) * m * (fin_span**2 + c_avg**2)
    Iyy = (1/12) * m * (fin_span**2 + fin_thickness**2)
    Izz = (1/12) * m * (c_avg**2 + fin_thickness**2)

    # local (body-frame-aligned) inertia tensor
    I_local = np.diag([Ixx, Iyy, Izz])

    # ------------------------------
    # 2. Build the multiple fins through rotation
    # ------------------------------
    masses = []
    for k in range(n_fins):
        angular_spacing = 360/n_fins
        angle = np.radians(angular_spacing * k)

        # Position CG around fuselage circumference
        cg = np.array([
            0.0,                     # x-axis along fuselage (place at x=0; you will shift later)
            0.5* od_fuse * np.cos(angle),  # y-axis
            0.5* od_fuse * np.sin(angle),  # z-axis
        ])

        # Rotation matrix to orient fin about x-axis
        R = np.array([
            [1,              0,              0],
            [0,  np.cos(angle), -np.sin(angle)],
            [0,  np.sin(angle),  np.cos(angle)]
        ])

        # Rotate inertia tensor into rocket coordinates
        I_rot = R @ I_local @ R.T

        # Create mass object
        masses.append(mass(mass=m,cg=cg,inertia=I_rot))

    return parallel_axis(masses) #returns mass obj for all fins





def tensor_to_rocketpy_inertia(tensor: np.ndarray):
    # convert numpy 3x3 matrix 
    return (float(tensor[0,0]), float(tensor[1,1]), float(tensor[2,2]))


def parallel_axis(mass_objs: list):
    """
    Input: list of mass objects

    Returns a mass object representing the total mass and inertia of the list 
    """

    mass_total = sum(mass_obj.mass for mass_obj in mass_objs)

    cg_total = sum(mass_obj.mass * np.array(mass_obj.cg) for mass_obj in mass_objs) / mass_total

    I_total = np.zeros((3, 3))
    for mass_obj in mass_objs:
        m = mass_obj.mass
        d = mass_obj.cg - cg_total
        d2 = np.dot(d, d)

        # Parallel axis contribution
        I_shift = mass_obj.inertia + m * (d2 * np.eye(3) - np.outer(d, d))
        I_total += I_shift

    return mass(mass_total, cg_total, I_total)



def empirical_cc_mass_model(V_cc, CC_LD, rho=4025):
    # Assume CC rho = 4025 kg/m^3, L/D ratio = 1.75 if none given
    
    d_outer = diameter_from_volume_LD(V_cc, CC_LD)

    L = CC_LD*d_outer

    print("d_outer: ", d_outer, L)

    d_inner = 0.8*d_outer #assume

    m = hollow_cylinder_mass(rho, L, d_inner, d_outer)

    return mass(m, np.array([(0.5*L),0.0,0.0]), hollow_cylinder_inertia(m, L, d_inner, d_outer) )




def empirical_rkt_mass_model(id_tank, od_tank, id_fuse, od_fuse, L_ox_tank, ox_tank_pos, L_nose, nose_position, rho_al, L_fuel_tank, fuel_tank_pos, rho_upperfuse, rho_lowerfuse, rho_nose, rho_fin, m_mev, m_ftv, m_otv, m_reco, fin_span, fin_root_chord, fin_tip_chord, fin_thickness, n_fins, h_nosecone):
    #PARSE THIS FUNCTION FOR ASSUMPTIONS, INPUT IN SI ONLY
    

    ### take point masses

    ### valve names are +x dir above the engine/tank, include tank bulkheads
    mev_pos = ox_tank_pos-0.5*L_ox_tank
    mev = mass(m_mev, np.array([mev_pos,0.0,0.0]), np.zeros((3,3)) )

    if L_fuel_tank != None:
        ftv_pos = fuel_tank_pos+L_fuel_tank
        ftv = mass(m_ftv, np.array([ftv_pos,0.0,0.0]), np.zeros((3,3)) )


    otv_pos = ox_tank_pos + 0.5*L_ox_tank
    otv = mass(m_otv, np.array([otv_pos,0.0,0.0]), np.zeros((3,3)) )


    #NOTE: reco might be a bit big to consider a pt mass

    if L_fuel_tank != None:
        reco_pos = (fuel_tank_pos + 0.5*L_fuel_tank) + 0.5*(nose_position - (fuel_tank_pos+L_fuel_tank) )
    else:
        reco_pos =  (ox_tank_pos + 0.5*L_ox_tank) + 0.5*(nose_position - (ox_tank_pos+L_ox_tank) )
    reco = mass(m_reco, np.array([reco_pos,0.0,0.0]), np.zeros((3,3)) )

    ### Cylinder Inertias
    ox_tank_masses = hollow_cylinder_lumped_masses(rho_al, L_ox_tank, id_tank, od_tank, ox_tank_pos)
    if L_fuel_tank != None:
        fuel_tank_masses = hollow_cylinder_lumped_masses(rho_al, L_ox_tank, id_tank, od_tank, fuel_tank_pos)



    if L_fuel_tank != None:
        L_upperfuse = (nose_position) - (0.5*L_fuel_tank + fuel_tank_pos)
        upperfuse_pos = 0.5*L_upperfuse + (0.5*L_fuel_tank + fuel_tank_pos)
    else:
        L_upperfuse = (nose_position) - (0.5*L_ox_tank + ox_tank_pos)
        upperfuse_pos = 0.5*L_upperfuse + (0.5*L_ox_tank + ox_tank_pos)
    upperfuse_masses = hollow_cylinder_lumped_masses(rho_upperfuse, L_upperfuse, id_tank, od_tank, upperfuse_pos)

    L_lowerfuse = (ox_tank_pos - 0.5*L_ox_tank)
    lowerfuse_pos = 0.5*L_lowerfuse
    lowerfuse_masses = hollow_cylinder_lumped_masses(rho_lowerfuse, L_lowerfuse, id_tank, od_tank, lowerfuse_pos)

    # fins and nosecone
    fins = total_fin_inertia(fin_span, fin_root_chord, fin_tip_chord, fin_thickness, rho_fin, od_fuse, n_fins)

    m_nosecone = thin_cone_shell_mass(rho_nose, h_nosecone, od_fuse)
    I_nosecone = thin_hollow_cone_inertia(m_nosecone, h_nosecone, od_fuse)
    #NOTE: HOW TO SOL NOSECONE POS SO ITS AT CENTROID, its just 1/3rd height + L to edge right?
    nosecone = mass(m_nosecone, np.array([nose_position,0.0,0.0]), I_nosecone)


    mass_data = [mev, otv, reco, fins, nosecone]
    for m in upperfuse_masses:
        mass_data.append(m)
    for m in lowerfuse_masses:
        mass_data.append(m)
    for m in ox_tank_masses:
        mass_data.append(m)

    if L_fuel_tank != None:
        mass_data.append(ftv)
        for m in fuel_tank_masses:
            mass_data.append(m)

    return mass_data






def mass_v1(id_tank, od_tank, id_fuse, od_fuse, L_ox_tank, ox_tank_pos, L_nose, nose_position, 
            rho_al, V_cc, rho_upperfuse, rho_lowerfuse, rho_nose, rho_fin, m_mev, m_ftv, m_otv, 
            m_reco, fin_span, fin_root_chord, fin_tip_chord, fin_thickness, n_fins, h_nosecone, 
            L_fuel_tank = None, fuel_tank_pos = None, CC_LD = 1.75):

    cc = empirical_cc_mass_model(V_cc, CC_LD)

    rkt_mass_data = empirical_rkt_mass_model(id_tank, od_tank, id_fuse, od_fuse, L_ox_tank, ox_tank_pos, L_nose, nose_position, rho_al, L_fuel_tank, fuel_tank_pos, rho_upperfuse, rho_lowerfuse, rho_nose, rho_fin, m_mev, m_ftv, m_otv, m_reco, fin_span, fin_root_chord, fin_tip_chord, fin_thickness, n_fins, h_nosecone)
    
    #plot rkt mass dist
    total_mass_data = rkt_mass_data.copy()
    total_mass_data.append(cc)
    plot_mass_distribution(total_mass_data)


    rkt = parallel_axis(rkt_mass_data)

    rkt_m = rkt.mass
    rkt_cg = rkt.cg[0]
    rkt_I = tensor_to_rocketpy_inertia(rkt.inertia)

    cc_m = cc.mass
    cc_cg = cc.cg[0]
    cc_I = tensor_to_rocketpy_inertia(cc.inertia)


    return rkt_m, rkt_cg, rkt_I, cc_m, cc_cg, cc_I, total_mass_data







def plot_mass_distribution(mass_objs, title="Mass Distribution"):
    """
    mass_objs: list of `mass` dataclass objects
    """

    xs = []
    ms = []
    labels = []

    for obj in mass_objs:
        xs.append(obj.cg[0])    # x-location
        ms.append(obj.mass)     # mass
        labels.append(str(obj)) # dataclass prints type info

    xs = np.array(xs)
    ms = np.array(ms)

    plt.figure(figsize=(10,4))
    sc = plt.scatter(xs, np.zeros_like(xs),
                     c=ms,
                     s=50 + 300*(ms/np.max(ms)),  # bubble size
                     cmap="viridis",
                     alpha=0.9)

    plt.axhline(0, color='black', linewidth=0.4)
    plt.xlabel("Position along rocket [m]")
    plt.title(title)

    # add colorbar for mass magnitude
    cb = plt.colorbar(sc)
    cb.set_label("Mass [kg]")

    # optional: annotate
    for x, m in zip(xs, ms):
        plt.text(x, 0.02, f"{m:.1f} kg", ha='center', fontsize=8)

    plt.tight_layout()
    plt.show()
