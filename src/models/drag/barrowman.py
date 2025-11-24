import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from dataclasses import dataclass
import math
from typing import Tuple

gamma = 1.4 #ROCKETPY USES CONSTANT RATIO OF SPECIFIC HEATS = 1.4

def sol_Cf_turbulent(Re, Rs, L_ref):
    
    Re_crit = 51*(Rs/L_ref)**(-1.039)

    if Re <= Re_crit:
        Cf = 1/((3.46*np.log10(Re)-5.6)**2) - 1700/Re
    else:
        Cf = 0.032*(Rs/L_ref)**2 #skin friction independent of Re

    return Cf

def load_C_DP_lookup(csv_path):
    df = pd.read_csv(csv_path)
    
    # Ensure sorted by Mach
    df = df.sort_values("Ma")

    x = df["Ma"].to_numpy()
    y = df["C_DP"].to_numpy()

    # Build linear interpolator with linear extrapolation
    f = interp1d(
        x,
        y,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate"
    )
    
    return f

def _busemann_coeffs(M: float, gamma: float=1.4) -> Tuple[float,float,float,float]:
    """
    Busemann second-order coefficients K1,K2,K3 and bow-shock term factor K* (per Eqs. (1)-(3)).
    K* is zero for purely expansive upstream flow; when a single upstream compressive shock exists,
    use Eq. (3). We return K1,K2,K3,B_shock where B_shock multiplies K*.
    """

    #NOTE: busemann doesnt work in transonic, and yea that is something this implementation does poorly but it will tie us over until we finalize config
    if M < 1.2:
        M = 1.2

    MSQ = M*M
    beta = math.sqrt(MSQ - 1.0)                           # Prandtl–Glauert factor
    # First- and second-order coefficients (Eq. (1) with (2a)-(2b))
    K1 = 2.0 / beta
    # K2 from Busemann theory with γ retained symbolically:
    K2 = (1.2*MSQ**2-2.0*beta**2)/beta**4 #((gamma+1.0)*M**4 - 4.0*beta**2) / (4.0*beta**4) 
    # Third-order coefficient K3 (Eq. (2c)). Barrowman keeps this term.
    DEM = -2.0*MSQ+4.0/3.0
    DEN = 1.0/beta**7
    K3 = (0.4*MSQ**4-(10.88/6.0)*MSQ**3+4.0*MSQ*MSQ+DEM)*DEN #(  (gamma+1.0)*M**8 + (2.0*gamma**2 -7.0*gamma -5.0)*M**6 + 10.0*(gamma+1.0)*M**4 -12*M**2 + 8.0 ) / (6*beta**7)

    # Bow-shock third-order irreversible contribution factor (multiplier of K*) – Eq. (3)
    Bshock = (0.04*MSQ**4-0.32*MSQ**3+0.4*MSQ*MSQ)*DEN#((gamma+1.0)*M**4 * ((5.0-3.0*gamma)*M**4 + 4.0*(gamma-3.0)*M**2 + 8.0)) / 48.0
    #print("Bshock: ", M, beta,K1,K2,K3,Bshock)
    return K1,K2,K3,Bshock

def _cp_local(eta: float, M: float, gamma: float,
              has_upstream_compressive_shock: bool) -> float:
    """
    Local Cp on a small flat region i at inclination η_i to the freestream,
    using Eq. (1): Cp = K1*η + K2*η^2 + K3*η^3 + K* (bow-shock term when applicable).
    Signs: positive η is a compressive (into the flow) tilt.
    """
    K1,K2,K3,Bshock = _busemann_coeffs(M, gamma)
    cp = K1*eta + K2*(eta**2) + K3*(eta**3)
    if has_upstream_compressive_shock:
        cp += Bshock * (eta**3)  # K* term – Eq. (1),(3) #NOTE: check this!
    return cp

def fin_aero_vs_alpha(
    M,
    a,
    A_ref,
    L_ref,
    C_root,
    S,
    Lambda_L,
    Lambda_1,
    Lambda_2,
    Lambda_T,
    zeta_L,
    zeta_T,
    lL_root,
    lT_root,
    n_strips,
    gamma,
):
    """
    Match FIN's primary output: for a given M, return CDw.
    All coefficients follow Barrowman's normalization for *two fins* for CL,CDw,CM,CLM;
    CPx (about reference axis) and CPS (about root chord) are for one fin (as in FIN).

    Geometry & marching:
      - Split span S into n strips of width Δy = S/n (Eq. (29)).
      - For each strip: build chord C(y) and region chords l_L(y), l_T(y), center region l_M (Eqs. (25)-(28),(30)-(32)).
      - Form the six local surface regions and their η_i (Eq. (33)), lengths d_i, n_i (Eq. (34)).
      - Apply Mach-cone correction as ratios r_i and ½·Cp behind the cone (Eqs. (36)–(44)).
      - Local lift/drag & x_cp via Eqs. (4)-(7); strip sums via (8)-(11); fin sums via (12)-(15).
      - Convert sums to coefficients via (16)-(19); compute CPx, CPS via (20)-(21).
    """

    # Precompute sweep tangents
    tL  = np.tan(np.deg2rad(Lambda_L))
    t1  = np.tan(np.deg2rad(Lambda_1))
    t2  = np.tan(np.deg2rad(Lambda_2))
    tT  = np.tan(np.deg2rad(Lambda_T))

    dy = S / float(n_strips)

    # Fin Tip Mach Cone Correction (Eq. (36))
    mu = math.asin(1.0/M)  # = arctan(1/beta) equivalently (by definition)

    a = np.deg2rad(a)
    TDRAG = 0.0

    y = 0.0
    for _ in range(n_strips):
        # Planform (Eqs. (25)-(28))
        XL = y * tL * math.cos(a)
        C  = y * (t2 - tL) + C_root
        lL = y * (t1 - tL) + lL_root
        lT = y * (tT - t2) + lT_root
        lM = max(C - lL - lT, 0.0)

        # Region geometry (Fig.4): 6 regions: (1,2,3) upper; (4,5,6) lower
        # Local surface inclinations η_i (Eq. (33a–f))
        etas = [
            #v tmp #NOTE
            np.deg2rad(zeta_L) - a,         # η1 = ζL − α
            -a,                           # η2 = −α
            -np.deg2rad(zeta_T) - a,        # η3 = −ζT − α
            np.deg2rad(zeta_L) + a,         # η4 = ζL + α
            +a,                           # η5 = +α
            -np.deg2rad(zeta_T) + a         # η6 = −ζT + α
        ]

        #for i in range(6):
        #    print("etas:", etas[i])
        #print("a, zeta: ", a, zeta_L, np.deg2rad(zeta_L))
        
        
        
        # Region surface lengths (Eqs. (30)–(32))
        Ls = [ lL, lM, lT, lL, lM, lT ]
        # Streamwise & normal components (Eq. (34))
        di = [ Ls[i]*math.cos(etas[i]) for i in range(6) ]
        ni = [ Ls[i]*math.sin(etas[i]) for i in range(6) ]
        
        #etas[0]*=-1.0

        #for i in range(6):
        #    print("direction:", di[i], ni[i], math.sqrt(di[i]**2 + ni[i]**2))
        
        # Forward region-boundary x-positions (Eq. (35))
        xp = [0.0]*6
        xp[0] = XL
        xp[1] = XL + di[0]
        xp[2] = XL + di[0] + di[1]
        xp[3] = XL
        xp[4] = XL + di[3]
        xp[5] = XL + di[3] + di[4]

        # Mach-cone correction ratios r_i (Eqs. (36)–(44))
        # distance from LE to cone cut on this strip (Eq. (38))
        #if include_mach_cone:

        lw = (S - y)*tL + np.tan(np.deg2rad(90-mu))
        ratios = [0.0]*6
        if lw >= C or lw <= 0.0:
            # no intersection or cone apex downstream of strip
            ratios = [0.0]*6
        elif lw >= C - lT:
            # only trailing region clipped (Eqs. (39)-(40))
            r_trail = (C - lw)/max(lT, 1e-16)
            ratios = [0.0, 0.0, 0.0, 0.0, 0.0, r_trail]
        elif lw >= lL:
            # center region clipped (Eqs. (41)-(42))
            r_mid = (C - lT - lw)/max(lM, 1e-16)
            ratios = [0.0, r_mid, 1.0, 0.0, r_mid, 1.0]
        else:
            # leading region clipped (Eq. (43)-(44))
            r_le = (lL - lw)/max(lL, 1e-16)
            ratios = [r_le, 1.0, 1.0, r_le, 1.0, 1.0]
        

        # Local Cp per region (Eq. (1)); flag K* term only for compressive side (upstream shock)
        # Barrowman’s program applies K* on the leading-edge surfaces (i=1 and i=4)

        cps = []
        for i in range(6):
            leading_face = (i in (0,3))
            cp = _cp_local(eta=etas[i], M=M, gamma=gamma, has_upstream_compressive_shock=leading_face)
            # Apply 1/2 coefficient in the area *behind* the Mach cone (Eq. (4)–(6) text)
            # Implement by scaling contribution, not Cp itself: (1 - 0.5*r_i)
            #print("cp: ", cp)
            cps.append(cp)


        # Local forces & moments on this strip (Eqs. (4)–(7); r_i handle behind-cone portion)
        FD_loc = 0.0
        for i in range(6):
            scale = (1.0 - 0.5*ratios[i])   # ½ Cp behind the cone
            FD_loc += cps[i] * ni[i] * scale


        # Accumulate strip → fin
        TDRAG += FD_loc

        y += dy

    # Convert to coefficients for two fins (Eqs. (16)–(19))
    # FIN divides by (q * A_ref) implicitly through Cp; here Cp already nondimensionalizes by q,
    # so the strip summations produce “per unit width”; multiply by Δy and normalize like FIN.
    # FIN’s implementation effectively multiplies by (Δy) and by number of fins as (2) or (4).
    # We keep the same scalings so the numerical outputs match the report.
    #CL  = 2.0 * TLIFT * dy / A_ref
    CDw = 2.0 * TDRAG * dy / A_ref

    #print("CDw: ", CDw, TDRAG, A_ref)

    return CDw



def make_fin_aero_solver(A_T,
            fuse_od,
            cr,
            S,
            Lambda_L,
            Lambda_1,
            Lambda_2,
            Lambda_T,
            zeta_L,
            zeta_T,
            lL_root,
            lT_root,
            n_strips: int = 300,
            gamma: float = 1.4):
    
    """
    Factory that returns a solver function depending ONLY on (M, alphas).
    All geometric & numerical settings are captured in the closure.
    """
    def solver(M: float, a):
        return fin_aero_vs_alpha(
            M=M,
            a=a,
            A_ref=A_T,
            L_ref=fuse_od,
            C_root=cr,
            S=S,
            Lambda_L=Lambda_L,
            Lambda_1=Lambda_1,
            Lambda_2=Lambda_2,
            Lambda_T=Lambda_T,
            zeta_L=zeta_L,
            zeta_T=zeta_T,
            lL_root=lL_root,
            lT_root=lT_root,
            n_strips=n_strips,
            gamma=gamma
        )
    return solver






def make_drag_func(fuse_od, nose_length, nose_position, fins_n, fins_span, fins_root_chord, fins_tip_chord, gamma_LE_sweep, tr, Lambda_L, Lambda_1, Lambda_2, Lambda_T, zeta_L, zeta_T, lL_root, lT_root):

    fuselage_radius = fuse_od/2

    #Constants: TODO: solve in factory function!
    L_ref =  nose_length+nose_position # ref length TODO: check tip to tail or body tube
    A_ref = 0.25*np.pi*fuse_od**2
    A_base = A_ref #NOTE: NO BOAT TAIL!
    Rs = 200e-6 #microns, incorrectly sprayed aircraft paint from barrowman thesis #20e-6 is well sprayed


    ### Fin Geometry
    #tr from inputs, root fin thickness
    #gamma_LE_sweep from inputs, midchord sweep angle
    S = fins_span                                   # Exposed fin semispan meas from root chord (barrowman and rocketpy naming does not agree, checked this)
    cr = fins_root_chord                            # fin root chord
    hr = tr*0.05                                    # fin root trailing edge thickness #NOTE: 0.05 is a rule of thumb for trailing edge fillet
    gamma_LE_sweep = np.deg2rad(gamma_LE_sweep)
    gamma_midchord = np.arctan( ( S*np.tan(gamma_LE_sweep) + (fins_tip_chord-cr)/2 )/S)  # leading edge sweep angle
    N = fins_n          # num fins
    A_Bf = 0.9*tr*cr             # base area of one fin #NOTE: 0.9 accounts for LE and TE fillet

    A_T = (cr+fins_tip_chord)/2*S                       # Planform area of one exposed fin panel
    rL =5e-4                                            # fin leading edge radius
    hf=2*rL                                             # Fin trailing edge thickness #NOTE: assume it has the same edge length as the fillet described by rL
    A_WB = (2*np.pi*fuselage_radius)*nose_position      # Wetted body area
    AR = (S**2)/A_T                                     # fin aspect ratio 


    
    ### Nose / Fuselage
    N_L = nose_length
    N_FR = nose_length/(2*fuselage_radius)                   # Nose finesse ratio
    FR = (nose_length+nose_position)/(2*fuselage_radius)     # Finesse ratio
    nosecone_area = 2*np.pi*fuselage_radius*nose_length*1.1
    A_w = N*A_T + A_WB + nosecone_area                       # Total Wetted Area

    csv_path = "src/models/drag/data/C_DP_TN_R100_Von_Karman.csv"
    C_DP_lookup = load_C_DP_lookup(csv_path)

    C_DWT_sol = make_fin_aero_solver(
        A_T,
        fuse_od,
        cr,
        S,
        Lambda_L,
        Lambda_1,
        Lambda_2,
        Lambda_T,
        zeta_L,
        zeta_T,
        lL_root,
        lT_root,

    )
    #print("C_DWT_sol inputs:",A_T,fuse_od,cr,S,Lambda_L,Lambda_1,Lambda_2,Lambda_T,zeta_L,zeta_T,lL_root, lT_root, )

    # The callable func:
    def drag_func(Ma, U, mu, rho, pressure, alpha, beta):

        if U < 0.01:
            U = 0.01 #Perturb for starting inst

        Re = rho*U*L_ref/mu
        Re_trans = 5e5


        eta = np.cbrt(AR*(tr/cr))
        C_DTT_crit = 1.15*((tr/cr)**(5/3))*(1.61+eta-np.sqrt( (eta-1.43)**2 + 0.578))

        C_DP_crit_ogive = 0.88/(N_FR+1)**2.22

        if Ma < 1.0: #Subsonic

            if Re < Re_trans: #Laminar

                Cf = 1.328/np.sqrt(Re)
                Cfc = Cf*(1-0.09*Ma**2)

                if Re < Re_trans: #Laminar
                    Cfc_crit = Cf*(0.91)
                else: # Re >= Re_trans, Turbulent
                    Cfc_crit = Cf*(0.88)

            else: # Re >= Re_trans, Turbulent

                Cf = sol_Cf_turbulent(Re, Rs, L_ref)
                Cfc = Cf*(1-0.12*Ma**2)

                if Re < Re_trans: #Laminar
                    Cfc_crit = Cf*(0.91)
                else: # Re >= Re_trans, Turbulent
                    Cfc_crit = Cf*(0.88)

            if Ma < 0.9:
                delta_CD = (1-Ma**2)**-0.417

            else: # 0.9 <= Ma < 1.0
                delta_CD = 1-1.5*(Ma-0.9)

            CfB = 2*Cfc*(cr/hr)
            K_DBT = (np.cos(np.deg2rad(gamma_midchord))**2) + ((0.223+4.02*Cfc_crit*(tr/hr)**2)**2)/np.cbrt((cr/hr)*Cfc_crit)**2

            #eqn 4-28
            C_DBT = 0.135*N*(A_Bf/A_ref)/( np.cbrt(CfB)*np.sqrt(K_DBT-(Ma**2)*np.cos(gamma_midchord)) )
                    
            K_DTT = (np.cos(np.deg2rad(gamma_midchord))**2) * np.cbrt( (C_DTT_crit/Cfc_crit)*(A_ref/A_T) -4*(tr/cr)*np.cos(np.deg2rad(gamma_midchord))/120*((tr/cr)**4)*(np.cos(np.deg2rad(gamma_midchord))**2) )**2

            #eqn 4-32
            C_DTT = 4*N*Cfc*(A_T/A_ref)*( (tr/cr)*np.cos(np.deg2rad(gamma_midchord)) + (30*((tr/cr)**4)*(np.cos(np.deg2rad(gamma_midchord))**2))/np.sqrt(K_DTT-(Ma**2)*(np.cos(np.deg2rad(gamma_midchord))**2))**3 )
            Cd_tail = C_DTT

            # Subsonic Body Drag Contribution

            K_DP = 1 + (6*A_WB*Cfc_crit)/( (FR**3)*A_ref*C_DP_crit_ogive)
            #eqn 4-41
            C_DP = (6*A_WB*Cfc)/( (FR**3)*A_ref*(K_DP-Ma**2)**0.6)

            K_DBB = 1 + 1/( (6.38 + 39.7*(hr/cr))*Cfc_crit*(A_w/A_base))
            #eqn 4-51
            C_DBB = (0.29*(A_base/A_ref))/np.sqrt(Cfc*(A_w/A_base)*(K_DBB-Ma**2))


    
        else: #Supersonic

            # Solve Laminar Case
            Cf = 1.328/np.sqrt(Re)
            #print("Re: ", Re)
            K = 0.15                        #K = 0.052 for cooled wall, K = 0.15 for no Q
            Cfc = Cf/( (1+K*Ma**2)**0.58)

            if Re >= Re_trans: #Turbulent

                Cf_turb = sol_Cf_turbulent(Re, Rs, L_ref)
                Cfc_turb = Cf_turb/(1+0.18*Ma**2)

                if Cfc_turb > Cfc:       # turbulent vals cant be less than laminar
                    Cf = Cf_turb
                    Cfc = Cfc_turb

            delta_CD = 1.214 - 0.502/Ma**2 + 0.1095/Ma**4 + 0.0231/Ma**6

            #eqn 4-30
            C_DBT = ( N*(1-0.52*Ma**-1.19)*(A_Bf/A_ref))/( (1+18*Cfc*(tr/hr)**2)*Ma**2)

            C_DWT = C_DWT_sol(Ma, alpha) + C_DWT_sol(Ma, beta) #(tr/cr) / np.sqrt( (Ma**2)*(np.cos(np.deg2rad(gamma_midchord))**2) -1) * (A_T/A_ref)

            #NOTE: SEEMS LIKE -1 IN SQRT IS WRONG
            #print("alpha and beta component: ", C_DWT_sol(Ma, alpha), C_DWT_sol(Ma, beta) )

            #print("C_DWT: ", C_DWT, alpha, beta, Ma )

            Cd_tail = C_DWT 

            # Supersonic Body Drag Contribution
            #print("check rktpy mod pressure: ", pressure)

            C_DP = C_DP_lookup(Ma)

            C_DB_crit = (A_base/A_ref)*(0.185+1.15*(hf/cr))
            Ma_crit = 0.892/np.sqrt(C_DB_crit)

            if Ma <= Ma_crit:
                C_DBB = C_DB_crit*( 0.88+0.12*np.exp(-3.58*(Ma-1) ))
            else:
                C_DBB = 0.7*(A_base/A_ref)/(Ma**2)


        #eqn: 4-15
        C_DFT = 2*N*Cfc*(A_T/A_ref)

        #eqn 4-22
        C_DLT = 2*N*(S*rL/A_ref)*(np.cos(np.deg2rad(gamma_LE_sweep))**2)*delta_CD


        Cd_tail += C_DFT + C_DLT + C_DBT
        Cd_body = C_DP + C_DBB

        Cd = Cd_tail + Cd_body

        #if Ma >= 1.0: #Subsonic
            #print("Cd!!", Cd, C_DWT, C_DFT + C_DLT + C_DBT, C_DP, C_DBB, Re)
        return Cd

    return drag_func
