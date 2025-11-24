import numpy as np
from dataclasses import dataclass

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

def thin_hollow_cone_mass(rho, H, d_inner, d_outer):
    
    Ri = d_inner / 2
    Ro = d_outer / 2

    # Thin shell area
    L = np.sqrt(H**2 + Ro**2)  # slant height (use outer R)
    A_shell = np.pi * (Ro + Ri) * L

    # Mass
    return rho * A_shell

def thin_hollow_cone_inertia(m, H, R):
    Ixx = 0.5 * m * R**2
    Iyy = 0.25 * m * R**2 + (1/12) * m * H**2
    Izz = Iyy
    inertia = np.array([[Ixx,0,0],[0,Iyy,0],[0,0,Izz]])
    return inertia



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




def empirical_rkt_mass_model(m_empirical_total, id_tank, od_tank, id_fuse, od_fuse, L_ox_tank, ox_tank_pos, L_nose, nose_position, rho_al, L_fuel_tank, fuel_tank_pos, m_mev = 13, m_ftv = 1.5, m_otv = 4.5, m_reco = 2.5):
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
    if L_fuel_tank != None:
        m_fuel_tank = hollow_cylinder_mass(rho_al, L_fuel_tank, id_tank, od_tank)
        I_fuel_tank = hollow_cylinder_inertia(m_fuel_tank, L_fuel_tank, id_tank, od_tank)
        fuel_tank = mass(m_fuel_tank, np.array([fuel_tank_pos,0.0,0.0]), I_fuel_tank)


    m_ox_tank = hollow_cylinder_mass(rho_al, L_ox_tank, id_tank, od_tank)
    I_ox_tank = hollow_cylinder_inertia(m_ox_tank, L_ox_tank, id_tank, od_tank)
    ox_tank = mass(m_ox_tank, np.array([ox_tank_pos,0.0,0.0]), I_ox_tank)

    ### sol remaining mass

    m_no_aerostruct = m_mev + m_otv + m_reco + m_ox_tank
    if L_fuel_tank != None:
        m_no_aerostruct += (m_fuel_tank + m_ftv)
    m_remaining = m_empirical_total - m_no_aerostruct

    #print("m_remaining:", m_remaining, m_empirical_total, - m_no_aerostruct)


    if L_fuel_tank != None:
        L_upperfuse = (0.5*L_nose + nose_position) - (L_fuel_tank + fuel_tank_pos)
        upperfuse_pos = 0.5*L_upperfuse + (L_fuel_tank + fuel_tank_pos)
    else:
        L_upperfuse = (0.5*L_nose + nose_position) - (L_ox_tank + ox_tank_pos)
        upperfuse_pos = 0.5*L_upperfuse + (L_ox_tank + ox_tank_pos)
    m_upperfuse = (2/5)*m_remaining #NOTE: ASSUMPTION, TODO: replace with physics based fin mass
    I_upperfuse = hollow_cylinder_inertia(m_upperfuse, L_upperfuse, id_fuse, od_fuse)
    upperfuse = mass(m_upperfuse, np.array([upperfuse_pos,0.0,0.0]), I_upperfuse)



    L_lowerfuse = (ox_tank_pos - 0.5*L_ox_tank) + 0
    lowerfuse_pos = 0.5*L_lowerfuse
    m_lowerfuse = (3/5)*m_remaining #NOTE: ASSUMPTION, TODO: replace with physics based fin mass
    I_lowerfuse = hollow_cylinder_inertia(m_lowerfuse, L_lowerfuse, id_fuse, od_fuse)
    lowerfuse = mass(m_lowerfuse, np.array([lowerfuse_pos,0.0,0.0]), I_lowerfuse)

    if L_fuel_tank != None:
        mass_data = [mev, ftv, otv, reco, fuel_tank, ox_tank, upperfuse, lowerfuse]
    else:
        mass_data = [mev, otv, reco, ox_tank, upperfuse, lowerfuse]

    #for m in mass_data:
    #    print(m)

    return mass_data



def mass_v1(m_empirical_total, id_tank, od_tank, id_fuse, od_fuse, L_ox_tank, ox_tank_pos, L_nose, nose_position, rho_al, V_cc, L_fuel_tank = None, fuel_tank_pos = None, CC_LD = 1.75):

    cc = empirical_cc_mass_model(V_cc, CC_LD)
    rktpy_m_empirical_total = m_empirical_total - cc.mass #dont double count cc!

    rkt_mass_data = empirical_rkt_mass_model(rktpy_m_empirical_total, id_tank, od_tank, id_fuse, od_fuse, L_ox_tank, ox_tank_pos, L_nose, nose_position, rho_al, L_fuel_tank, fuel_tank_pos)
    
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


    return rkt_m, rkt_cg, rkt_I, cc_m, cc_cg, cc_I





import matplotlib.pyplot as plt
import numpy as np

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
