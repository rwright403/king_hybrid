import numpy as np
gamma = 1.4 #ROCKETPY USES CONSTANT RATIO OF SPECIFIC HEATS = 1.4

def sol_Cf_turbulent(Re, Rs, L_ref):
    
    Re_crit = 51*(Rs/L_ref)**(-1.039)

    if Re <= Re_crit:
        Cf = 1/((3.46*np.log10(Re)-5.6)**2) - 1700/Re
    else:
        Cf = 0.032*(Rs/L_ref)**2 #skin friction independent of Re

    return Cf



def make_drag_func(fuse_od, nose_length, nose_position, fins_n, fins_span, fins_root_chord, fins_tip_chord, gamma_LE_sweep, tr):

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


    # The callable func:
    def drag_func(Ma, U, mu, rho, pressure):

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
            print("Re: ", Re)
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

            C_DWT = N * (tr/cr) / np.sqrt( (Ma**2)*(np.cos(np.deg2rad(gamma_midchord))**2) -1) * (A_T/A_ref)

            Cd_tail = C_DWT 

            # Supersonic Body Drag Contribution

            C_DP = cdp_supersonic_rkt(Ma, pressure, N_L, fuselage_radius)

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

        print("Cd!!", Cd, C_DP, C_DBB, Re)
        return Cd

    return drag_func

