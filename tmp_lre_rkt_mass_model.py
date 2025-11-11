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



def empirical_lre_cc_mass_model(V_cc, rho=4025, LD = 1.75):
    # Assume CC rho = 4025 kg/m^3, L/D ratio = 1.75 if none given
    
    d_outer = diameter_from_volume_LD(V_cc, LD)
    L = LD*d_outer

    d_inner = 0.6*d_outer #assume

    m = hollow_cylinder_mass(rho, L, d_inner, d_outer)
    tensor = hollow_cylinder_inertia(m, L, d_inner, d_outer)

    return tensor_to_rocketpy_inertia(tensor)


def empirical_rkt_mass_model(m_empirical_total, id_tank, od_tank, id_fuse, od_fuse, L_ox_tank, ox_tank_pos, L_nose, nose_position, rho_al, L_fuel_tank = None, fuel_tank_pos = None, m_mev = 4, m_ftv = 1.5, m_otv = 2.5, m_reco = 4.5):
    #PARSE THIS FUNCTION FOR ASSUMPTIONS, INPUT IN SI ONLY

    ### take point masses

    ### valve names are +x dir above the engine/tank, include tank bulkheads
    mev_pos = ox_tank_pos-0.5*L_ox_tank
    mev = mass(m_mev, np.ndarray([mev_pos,0.0,0.0]), np.zeros((3,3)) )

    ftv_pos = fuel_tank_pos+L_fuel_tank
    ftv = mass(m_ftv, np.ndarray([ftv_pos,0.0,0.0]), np.zeros((3,3)) )

    otv_pos = (ox_tank_pos + 0.5*L_ox_tank) + 0.5*( (fuel_tank_pos-0.5*L_fuel_tank) - (ox_tank_pos + 0.5*L_ox_tank) )
    otv = mass(m_otv, np.ndarray([otv_pos,0.0,0.0]), np.zeros((3,3)) )

    #NOTE: reco might be a bit big to consider a pt mass
    reco_pos = (fuel_tank_pos + 0.5*L_fuel_tank) + 0.5*(nose_position - (fuel_tank_pos+L_fuel_tank) )
    reco = mass(m_reco, np.ndarray([reco_pos,0.0,0.0]), np.zeros((3,3)) )

    ### Cylinder Inertias

    L_fuel_tank #passed in
    m_fuel_tank = hollow_cylinder_mass(rho_al, L_fuel_tank, id_tank, od_tank)
    I_fuel_tank = hollow_cylinder_inertia(m_fuel_tank, L_fuel_tank, id_tank, od_tank)
    fuel_tank = mass(m_fuel_tank, np.ndarray([fuel_tank_pos,0.0,0.0]), I_fuel_tank)


    m_ox_tank = hollow_cylinder_mass(rho_al, L_ox_tank, id_tank, od_tank)
    I_ox_tank = hollow_cylinder_inertia(m_ox_tank, L_ox_tank, id_tank, od_tank)
    ox_tank = mass(m_ox_tank, np.ndarray([ox_tank_pos,0.0,0.0]), I_ox_tank)

    ### sol remaining mass

    m_no_aerostruct = m_mev + m_otv + m_ftv + m_reco + m_fuel_tank + m_ox_tank
    m_remaining = m_empirical_total - m_no_aerostruct



    L_upperfuse = (0.5*L_nose + nose_position) - (L_fuel_tank + fuel_tank_pos)
    upperfuse_pos = 0.5*L_upperfuse + (L_fuel_tank + fuel_tank_pos)
    m_upperfuse = (2/5)*m_remaining #NOTE: ASSUMPTION, TODO: replace with physics based fin mass
    I_upperfuse = hollow_cylinder_inertia(m_upperfuse, L_upperfuse, id_fuse, od_fuse)
    upperfuse = mass(m_upperfuse, np.ndarray([upperfuse_pos,0.0,0.0]), I_upperfuse)

    L_lowerfuse = (ox_tank_pos - 0.5*L_ox_tank) + 0
    lowerfuse_pos = 0.5*L_lowerfuse
    m_lowerfuse = (3/5)*m_remaining #NOTE: ASSUMPTION, TODO: replace with physics based fin mass
    I_lowerfuse = hollow_cylinder_inertia(rho_al, L_lowerfuse, id_fuse, od_fuse)
    lowerfuse = mass(m_lowerfuse, np.ndarray([lowerfuse_pos,0.0,0.0]), I_lowerfuse)


    mass_data = [mev, ftv, otv, reco, fuel_tank, ox_tank, upperfuse, lowerfuse]
    dry_rkt_no_motor = parallel_axis(mass_data)

    return tensor_to_rocketpy_inertia(dry_rkt_no_motor.inertia)



