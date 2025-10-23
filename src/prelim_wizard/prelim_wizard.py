from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel
import CoolProp.CoolProp as CP #I love coolprop! ~ units: http://www.coolprop.org/v4/apidoc/CoolProp.html
import numpy as np
import matplotlib.pyplot as plt


def bar_to_psi(x):
    return x * 14.5038


def magic(inputs):
    """
    Initial Design Spreadsheet but in python

    See total impulse est It_est: 
    Used other team's rockets to create a relationship (linear regression) between:
    apogee, total impulse and wet mass
    (search for DEV LRE Hand Calcs on the team google drive)   
    """

    ### Step 0:
    #estimate vapor pressure wrt conditions
    if inputs.oxidizer_name  == "N2O":

        # Define temperature range (input in celsius, program converts to kelvin)
        T_min = 16  # Minimum operating temperature
        T_max = 27  # Maximum operating temperature
        T_crit = 36 # Critical temperature of nitrous
        num_points = 100  # Number of points for the plot

        # Generate an array of temperatures
        T_min = T_min + 273.15
        T_max = T_max + 273.15
        T_crit = T_crit + 273.15
        temperatures = np.linspace(T_min, T_crit, num_points)
        temperatures_celsius = temperatures - 273.15  # Convert from Kelvin to Celsius
        vapor_pressures = np.zeros_like(temperatures)

        # Plotting
        plt.figure(figsize=(10, 6))

        # Calculate vapor pressure for each temperature using CoolProp
        for i, T in enumerate(temperatures):
            P = CP.PropsSI('P', 'T', T, 'Q', 0, 'N2O')
            vapor_pressures[i] = P / 100000.0  # Convert from Pa to bar

        # Plot vapor pressure for the current quality
        plt.plot(temperatures_celsius, vapor_pressures, color='k', label='N2O Saturated Liquid Line')

        # THESE TEMPERATURES ARE MAX AND MIN VALUES FOR HOTFIRING
        plt.axvline(x=T_min-273.15, color='b', linestyle='--', label=f'{T_min-273.15} °C')
        plt.axvline(x=T_max-273.15, color='orange', linestyle='--', label=f'{T_max-273.15} °C')
        plt.axhline(y=72.5, color='r', linestyle='--', label="N2O Critical Pressure")

        plt.xlabel('Temperature (°C)')
        plt.ylabel('Vapor Pressure (bar)')
        plt.title('Vapor Pressure of Nitrous Oxide (N2O) vs. Temperature')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()


    ### Setup
    optimal_height = (2/3)*inputs.apogee_height #m above launch pad
    
    """
    #NOTE: this uses a rule of thumb, works rlly well for small change in alt sounding rocket flight (10k ft)
    does not work as well for 30k ft

    However, this is just an initial guess to approximate P_exit so we can solve for some things, and then resolve it when we have a better idea (circular calculations)
    """
    P_exit = 1000 * 101.29*( ((15.04 - 0.00649*optimal_height)+273.1)/288.08)**5.256 #Pa


    #NOTE: this is approx the range we should see for a student pressure fed nitrous oxide engine, change as required
    OF_ratio = [4, 5, 6, 7, 8, 9]
    chamber_pressures = [10, 20, 30, 40, 50]




    ### Step 1: setup Pcc, O/F, and combustion property plots to select P_cc and O/F

    # Create CEA object for the first graph (Flame Temp/mw vs. Chamber Pressure)
    ceaObj = CEA_Obj(
        oxName=inputs.oxidizer_name, fuelName=inputs.fuel_name, pressure_units='Bar', isp_units='sec', cstar_units='m/s',
        temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='kJ/kg', density_units='kg/m^3',
        specific_heat_units='kJ/kg-K'
    )


    # Create subplots
    fig, axs = plt.subplots(1, 3, figsize=(15, 6))

    # First subplot: Flame Temp/mw vs. Chamber Pressure
    axs[0].set_title('Flame Temp/mw vs Chamber Pressure')
    axs[0].set_xlabel('Chamber Pressure (Bar)')
    axs[0].set_ylabel('Flame Temp/mw ((mol K)/kg)')

    axs[1].set_title('Flame Temp vs Chamber Pressure')
    axs[1].set_xlabel('Chamber Pressure (Bar)')
    axs[1].set_ylabel('Flame Temp (K)')

    for j in OF_ratio:
        flame_temp_mw_ratio_arr = []
        flame_temp_arr = []
        p_cc_arr = []

        for k in chamber_pressures:

            ###solve expansion ratio for each chamber pressure using y guess
            ###pick ideal expansion altitude and solve expansion ratio

            #first need to get ratio of specific heats (y), from testing noticed in range of expected chamber pressures
            #and expansion ratios y was constant to expratio and only changed with P_cc so we can solve here w/out losing accuracy

            y = 1.23 #fluid_prop[1] # (-)

            expratio = ( ( ((y+1)/2)**(1/(y-1)) ) * ( (P_exit/k)**(1/y) ) * np.sqrt( ((y+1)/(y-1)) * ( (1- (P_exit/k)**((y-1)/y) )) ) )**-1


            i = ceaObj.get_IvacCstrTc_ChmMwGam(Pc=k, MR=j, eps=expratio)
            flame_temp_mw_ratio_arr.append(i[2] / i[3])
            flame_temp_arr.append(i[2])
            p_cc_arr.append(k)

        print(f"O/F Ratio: {j}, gamma = {i[4]}, gas const = {8134 / i[3]}, mw = {i[3]}")
        
        #graph OF vs flame temp / mw
        axs[0].plot(p_cc_arr, flame_temp_mw_ratio_arr, label=f'O/F={j}')
        #graph OF vs flame temp
        axs[1].plot(p_cc_arr, flame_temp_arr, label=f'O/F={j}')

    axs[0].legend()
    axs[0].grid(True)
    axs[1].legend()
    axs[1].grid(True)

    # Second subplot: ISP vs. O/F Ratio at Different Chamber Pressures


    axs[2].set_title('ISP vs O/F Ratio at Different Chamber Pressures')
    axs[2].set_xlabel('O/F Ratio')
    axs[2].set_ylabel('ISP (s)')

    C = CEA_Obj(
        oxName=inputs.oxidizer_name, fuelName=inputs.fuel_name, pressure_units='Pa', isp_units='sec', cstar_units='m/s',
        temperature_units='K', sonic_velocity_units='m/s', enthalpy_units='kJ/kg', density_units='kg/m^3',
        specific_heat_units='kJ/kg-K'
    )

    for Pc in chamber_pressures:


        of_arr = []
        isp_arr = []

        expratio = 5.5
        isp_prev = 0
        max_isp = 0
        of_store = 0
        i = 0.55

        while i < 16:
            isp = C.get_Isp(Pc=(Pc)*1e5, MR=i, eps=expratio, frozen=0, frozenAtThroat=0)
            of_arr.append(i)
            isp_arr.append(isp)

            if isp > isp_prev:
                max_isp = isp
                of_store = i
            i = i + 0.5
            isp_prev = isp

        print(f"For Pc={Pc} bar, Max ISP: {max_isp}, O/F Ratio: {of_store}")
        #graph OF vs max theoretical ISP
        axs[2].plot(of_arr, isp_arr, label=f'Pc={Pc} bar [{bar_to_psi(Pc)} PSI]')

    axs[2].legend()
    axs[2].grid(True)

    print("Pick O/F and P_cc, close graph and enter when prompted:")
    plt.show()

    selected_OF = 0
    while(selected_OF == 0):
        selected_OF = float(input("Enter O/F ratio: "))

    selected_Pcc = 0
    while(selected_Pcc == 0):
        selected_Pcc = float(input("Enter P_cc (bar): "))

    selected_Pcc *= 1e5         ###NOTE: UNIT CONVERSION FROM BAR TO Pa!!!!!






    ### Step 2: solving the throat area

    #use atmospheric model to get exit pressure at optimal height
    #NOTE: This is NASA Std. Atmospheric Model

    #P_pad = 1000 * 101.29*( ((15.04 - 0.00649*0)+273.1)/288.08)**5.256 #Pa #NOTE: sea level for getting base program
    P_exit = 1000 * 101.29*( ((15.04 - 0.00649*optimal_height)+273.1)/288.08)**5.256 #Pa
    #P_apogee = 1000 * 101.29*( ((15.04 - 0.00649*inputs.apogee_height)+273.1)/288.08)**5.256 #Pa #NOTE: sea level for getting base program


    #solve expansion ratio for optimal height
    #print(P_exit, selected_Pcc)
    expratio = ( ( ((y+1)/2)**(1/(y-1)) ) * ( (P_exit/selected_Pcc)**(1/y) ) * np.sqrt( ((y+1)/(y-1)) * ( (1- (P_exit/selected_Pcc)**((y-1)/y) )) ) )**-1
    #print("expansion ratio: ", expratio)


    #solve inital and final Cf, assuming rocket can always reach apogee (uses expansion ratio!!!!)
    Cf_opt = np.sqrt( ((2*y**2)/(y-1)) * ( (2/(y+1))**((y+1)/(y-1)) ) * (1- (P_exit/selected_Pcc)**((y-1)/y)) ) 

    #NOTE: start calculate impulse and rocket dry mass with spreadsheet line of best fit eqn
    
    if hasattr(inputs, "fuel_tank_model"):
        """
        LRE Linear Regression Estimates from spreadsheet
        """
        It_est = 2.73*inputs.apogee_height + 4829
        rocket_dry_mass = (1.03e-3)*It_est + 21
    else:
        """
        HYBRID Linear Regression Estimates from spreadsheet
        """
        It_est = 2.41*inputs.apogee_height + 11942
        rocket_dry_mass = (8.05e-5)*It_est + 54.7

    #use impulse estimation to graph a bunch of preliminary thrust curves that differ based on burn time and show the user
    burn_time_arr = [2, 3, 4, 5, 6, 7, 8, 9, 10]

    for t in burn_time_arr:

        #solve throat area
        A_throat = It_est / (t*selected_Pcc*Cf_opt)

        #solve thrust
        F_thrust = selected_Pcc*A_throat*Cf_opt

        #plot curve
        
        time = [0, t, t]  # Example times (t0, t1, t2)
        thrust = [F_thrust, F_thrust, 0]  # Example thrust values (initial, final, 0)

        # Plot the thrust curve
        plt.plot(time, thrust, label=f't={t}, Throat Diam={39.3701*(np.sqrt(4*A_throat/np.pi))} (in)')


    min_start_thrust = (rocket_dry_mass*9.81) * inputs.min_TW_ratio
    print(f"For reference: to achieve the estimated min_start_thrust of {inputs.min_TW_ratio}, need a starting thrust of {min_start_thrust:.2f} (N).")

    plt.title(label=f'Preliminary Thrust Curve Based on Estimated Total Impulse' )
    plt.xlabel('Burn Time (s)')
    plt.ylabel('Thrust (N)')
    plt.grid(True)
    plt.legend()
    plt.show()


    #allow the user to pick the burn time
    selected_tburn = 0
    while(selected_tburn == 0):
        selected_tburn = float(input("Pick and Enter Burn Time (s): "))






    ### Step 3: Evaluate Nozzle Performance
    A_throat = It_est / (selected_tburn*selected_Pcc*Cf_opt)
    throat_diam = 39.3701*(2 *np.sqrt(A_throat/np.pi))

    """ #OLD: USES RULE OF THUMB TO OPTIMIZE TO 2/3RDS EXIT AREA
    #solve exit area
    expratio = ( ( ((y+1)/2)**(1/(y-1)) ) * ( (P_exit/selected_Pcc)**(1/y) ) * np.sqrt( ((y+1)/(y-1)) * ( (1- (P_exit/selected_Pcc)**((y-1)/y) )) ) )**-1
    print("expansion ratio: ", expratio)
    A_exit = expratio * A_throat


    height_arr = np.linspace(inputs.elevation, (inputs.elevation+inputs.apogee_height), 50)
    pressure_arr = []
    for h in height_arr:
        t = 15.04 - 0.00649*h #Celsius
        p = 1000 * 101.29*( (t+273.1)/288.08)**5.256 #Pa
        pressure_arr.append(p) #Pa

    thrust_arr = []
    optimal_thrust_arr = []

    for P_atm in pressure_arr:

        #solve equation 1 (need to call rocketcea to get ratio of specific heats)
        cf = np.sqrt( ((2*y**2)/(y-1)) * ( (2/(y+1))**((y+1)/(y-1)) ) * (1- (P_exit/selected_Pcc)**((y-1)/y)) ) + (((P_exit-P_atm)*A_exit)/(selected_Pcc*A_throat))
        optimal_cf = np.sqrt( ((2*y**2)/(y-1)) * ( (2/(y+1))**((y+1)/(y-1)) ) * (1- (P_atm/selected_Pcc)**((y-1)/y)) )

        #solve equation 2 to get thrust
        thrust = cf * A_throat * selected_Pcc
        optimal_thrust = optimal_cf * A_throat * selected_Pcc

        #add thrust and convert exit pressure to altitude and add that to an array.
        thrust_arr.append(thrust)
        optimal_thrust_arr.append(optimal_thrust)

    #plot vertical line for minimum thrust to weight ratio
    min_start_thrust = (rocket_dry_mass*9.81) * inputs.min_TW_ratio
    plt.axvline(x=min_start_thrust, color='r', label=f'min starting thrust for T/W of {inputs.min_TW_ratio}')
    #plot

    print("Summary of Nozzle Performance:")
    plt.plot(thrust_arr, height_arr, label=f'optimal expansion at {optimal_height} (m)')
    plt.plot(optimal_thrust_arr, height_arr, label=f'optimal expansion throughout burn')
    plt.title(label=f'Altitude vs Thrust for throat_diam (in) {throat_diam}' )
    plt.xlabel('Thrust (N)')
    plt.ylabel('Altitude (m)')
    plt.grid(True)
    plt.legend()
    plt.show()
    """

    # === Step 3: Evaluate Nozzle Performance  (OPTIMIZED DESIGN ALTITUDE) ===
    # keep your existing A_throat sizing from It_est / selected_tburn / Cf_opt
    A_throat = It_est / (selected_tburn * selected_Pcc * Cf_opt)
    throat_diam = 39.3701 * (2.0 * np.sqrt(A_throat/np.pi))

    # Use your notation explicitly
    P_cc = float(selected_Pcc)      # chamber pressure [Pa]
    A_t  = float(A_throat)
    gamma = float(y)                # ratio of specific heats

    # Altitude grid for the area calculation
    h0 = float(inputs.elevation)
    H  = float(inputs.elevation + inputs.apogee_height)
    height_arr = np.linspace(h0, H, 300)

    # --- local helpers (fine to be nested) ---
    def P_atm(h):
        """ISA-ish model (Pa)."""
        t = 15.04 - 0.00649*h  # Celsius
        return 1000.0 * 101.29 * ((t + 273.1) / 288.08)**5.256

    def Cf_isent(P_exit_local):
        """Isentropic Cf for given exit static pressure (Pa)."""
        pr = np.clip(P_exit_local / P_cc, 1e-12, 0.999999)
        term = 1.0 - pr**((gamma - 1.0)/gamma)
        return np.sqrt((2*gamma**2/(gamma-1.0)) * (2/(gamma+1.0))**((gamma+1.0)/(gamma-1.0)) * term)

    def eps_from_Pexit(P_exit_local):
        """Expansion ratio eps = A_e/A_t for given P_exit and P_cc."""
        pr = np.clip(P_exit_local / P_cc, 1e-12, 0.999999)
        num = ((gamma+1.0)/2.0)**(1.0/(gamma-1.0)) * pr**(1.0/gamma)
        den = np.sqrt(((gamma+1.0)/(gamma-1.0)) * (1.0 - pr**((gamma-1.0)/gamma)))
        return 1.0 / (num * den)

    def build_curves(h_star):
        """
        Return thrust arrays for:
        - blue: fixed nozzle sized at design altitude h_star (P_exit_star = P_atm(h_star))
        - orange: always-optimally-expanded (P_exit = P_atm(h))
        Also returns eps(h_star).
        """
        P_exit_star = P_atm(h_star)         # design exit pressure (fixed nozzle)
        eps_star    = eps_from_Pexit(P_exit_star)
        A_e         = eps_star * A_t

        P_a = P_atm(height_arr)

        # Orange: always optimally expanded (P_exit = P_a)
        Cf_opt_all = Cf_isent(P_a)

        # Blue: single nozzle sized at h_star (add pressure mismatch term)
        Cf_fixed = Cf_isent(P_exit_star) + (P_exit_star - P_a)/P_cc * eps_star

        T_opt   = P_cc * A_t * Cf_opt_all
        T_fixed = P_cc * A_t * Cf_fixed
        return T_fixed, T_opt, eps_star, P_exit_star

    def objective(h_star, L2=False):
        """Area between curves over altitude (default L1)."""
        T_fixed, T_opt, _, _ = build_curves(h_star)
        diff = T_fixed - T_opt
        if L2:
            return np.trapz(diff**2, height_arr)
        else:
            return np.trapz(np.abs(diff), height_arr)

    # --- optimize design altitude h* to minimize area ---
    try:
        from scipy.optimize import minimize_scalar
        res = minimize_scalar(lambda h: objective(h, L2=False),
                            bounds=(h0, H), method="bounded",
                            options={"xatol": 1.0})
        h_star_opt = float(res.x)
        area_min   = float(res.fun)
    except Exception:
        # Fallback: grid search if SciPy isn't available
        candidates = np.linspace(h0, H, 401)
        vals = np.array([objective(h) for h in candidates])
        idx = int(np.argmin(vals))
        h_star_opt = float(candidates[idx])
        area_min   = float(vals[idx])

    # Build optimum curves and parameters
    T_fixed_opt, T_opt_opt, eps_opt, P_exit_star = build_curves(h_star_opt)
    A_exit = eps_opt * A_t
    expratio = eps_opt  # use optimized expansion ratio downstream

    # Plot vertical line for minimum T/W
    min_start_thrust = (rocket_dry_mass*9.81) * inputs.min_TW_ratio

    print("\n=== Nozzle design-altitude optimization ===")
    print(f"Optimal design altitude h*: {h_star_opt:.1f} m")
    print(f"Design exit pressure P_exit(h*): {P_exit_star:.0f} Pa")
    print(f"Expansion ratio eps = A_e/A_t: {eps_opt:.3f}")
    print(f"Exit area A_exit: {A_exit:.6f} m^2")
    print(f"Area between curves (L1): {area_min:.3e} N·m")
    print(f"Average mismatch per meter: {area_min/(H-h0):.2f} N")

    # Plot (blue = fixed @ h*, orange = always-optimal)
    plt.figure()
    plt.plot(T_fixed_opt, height_arr, label=f'fixed nozzle @ h*={h_star_opt:.0f} m')
    plt.plot(T_opt_opt,   height_arr, label='optimal expansion throughout burn')
    plt.axvline(x=min_start_thrust, color='r', label=f'min starting thrust for T/W of {inputs.min_TW_ratio}')
    plt.title(f'Altitude vs Thrust (throat_diam {throat_diam:.3f}\" ; eps={eps_opt:.2f})')
    plt.xlabel('Thrust (N)'); plt.ylabel('Altitude (m)')
    plt.grid(True); plt.legend(); plt.show()

    # For later summaries, use an average thrust from the optimized fixed-nozzle curve:
    thrust = float(np.trapz(T_fixed_opt, height_arr) / (height_arr[-1] - height_arr[0]))


    #solve exit velocity under those conditions
    fluid_prop = C.get_Chamber_MolWt_gamma(selected_Pcc, selected_OF, expratio)
    R = 8314 / fluid_prop[0] # J/(kg K)
    y = fluid_prop[1] # (-)

    temperatures = C.get_Temperatures(selected_Pcc, selected_OF, expratio, 0, 1)
    T_cc = temperatures[0]

    v_exit = np.sqrt(((2 * y) / (y - 1)) * R * T_cc * (1 - (P_exit / selected_Pcc)**((y - 1) / y)))

    m_dot_cc = thrust / v_exit

    m_dot_fuel = m_dot_cc/(selected_OF +1)
    m_dot_ox = m_dot_cc - m_dot_fuel

    m_propellant = rocket_dry_mass * inputs.mass_fraction_estimate

    m_fuel = m_propellant/(selected_OF +1)
    m_ox = m_propellant - m_fuel


    #pick tank pressures:
    selected_P_ox_tank = 0
    while(selected_P_ox_tank == 0):
        selected_P_ox_tank = float(input("Pick and Enter Oxidizer Tank Pressure (Bar): "))*1e5

    selected_P_fuel_tank = 0
    if hasattr(inputs, "fuel_tank_model"):
        while(selected_P_fuel_tank == 0):
            selected_P_fuel_tank = float(input("Pick and Enter Fuel Tank Pressure (Bar): "))*1e5


    #estimate injector area!!!!
    print(f"estimating reqiured injector areas for fuel and oxidizer assuming a Cd of {inputs.Cd_est} and SPI model")
    
    #estimate propellant densities with coolprop
    rho_oxidizer = CP.PropsSI('D', 'P', selected_P_ox_tank, 'T', 275, inputs.oxidizer_name)
    rho_fuel = 900 #kg/m
    if hasattr(inputs, "fuel_tank_model"):
        rho_fuel = CP.PropsSI('D', 'P', selected_P_fuel_tank, 'T', 275, inputs.fuel_name)
    
    A_ox_inj = m_dot_ox / (inputs.Cd_est*np.sqrt(2*rho_oxidizer*(selected_P_ox_tank-selected_Pcc)))
    A_fuel_inj = None
    if hasattr(inputs, "fuel_tank_model"):
        A_fuel_inj = m_dot_fuel / (inputs.Cd_est*np.sqrt(2*rho_fuel*(selected_P_fuel_tank-selected_Pcc)))
    
    ### Finally, print a summary of the results

    isp = C.get_Isp(Pc=selected_Pcc, MR=selected_OF, eps=expratio, frozen=0, frozenAtThroat=0)
    c_star = C.get_Cstar(Pc=selected_Pcc, MR=selected_OF)

    print(f"\n\n\n############\nPreliminary Design Summary:\n############\n\n")
    print(f"Performance:\n------------\nSpecific Impulse: {isp} (s)\nCharacteristic Velocity: {c_star} (m/s)\nBurn Time: {selected_tburn} (s)\nRocket Dry Mass: {rocket_dry_mass} (kg)\nTotal Impulse: {It_est} (N s)\n")
    print(f"Combustion Chamber:\n------------\nP_cc: {selected_Pcc} (Pa)\nO/F Ratio: {selected_OF}")
    print(f"Flame Temperature: {T_cc}\nRatio of Specific Heats: {y}\nReactant Gas Const. {R} (J/(kg K))")
    print(f"Average Thrust: {thrust} (N)\nMass Flow Rate: {m_dot_cc} (kg/s)\nExit Velocity: {v_exit} (m/s)")
    print(f"Expansion Ratio: {expratio}\nThroat Diameter {throat_diam} (in)\n")
    print(f'Feed System:\n------------\nEstimated Mass Fraction: {inputs.mass_fraction_estimate}\nFuel Mass: {m_fuel} (kg)\nOxidizer Mass: {m_ox} (kg)')
    print(f"Oxidizer Tank Pressure: {selected_P_ox_tank} (Pa)\nFuel Tank Pressure: {selected_P_fuel_tank} (Pa)")
    print(f"Preliminary Injector Solved by Assuming Fuel and Oxidizer Orifices Follow SPI Model and a Cd guess of {inputs.Cd_est}")
    print(f"Oxidizer Injector Discharge Area: {A_ox_inj} (m^2), Fuel Injector Discharge Area {A_fuel_inj} (m^2)\n")


    # forgot to size tanks my bad:
    if hasattr(inputs, "fuel_tank_model"):
        V_fuel = m_fuel / rho_fuel  # m^3

        #NOTE: BAD CODE INPUTS HERE AGILE
        ullage_fraction = 0.15
        fluid_gas = "Nitrogen"  # ullage gas

        # --- Compute total tank volume including ullage ---
        R_univ = 8.314462618  # J/(mol·K)
        M_molar = CP.PropsSI('M','Nitrogen')  # kg/mol
        R_specific = R_univ / M_molar      # J/(kg·K)

        V_tank = V_fuel / (1 - ullage_fraction)
        V_ullage = V_tank - V_fuel

        # --- Compute mass of pressurant using ideal gas law ---
        m_pressurant = selected_P_fuel_tank * V_ullage / (R_specific * 0.5*(T_min+T_max))

        # --- Temperature range ---
        T_range = np.linspace(273, 323, 100)  # K, e.g., 0°C to 50°C

        # --- Compute pressure at each temperature ---
        P = m_pressurant * R_specific * T_range / V_ullage  # Pa

        # --- Plot ---
        plt.figure(figsize=(8,5))
        plt.plot(T_range-273.15, P/1e5, color='blue')  # x in °C, y in bar
        plt.xlabel("Temperature [°C]")
        plt.ylabel("Nitrogen Pressure [bar]")
        plt.title("Ullage Gas Pressure vs Temperature")
        plt.grid(True)
        plt.show()




    if inputs.oxidizer_name  == "N2O":

        # Parameters
        fluid = "N2O"
        V_min = m_ox / CP.PropsSI("D","T",T_max,"Q",0,fluid)  # min tank volume
        V_max = V_min * 1.75
        V_points = np.linspace(V_min, V_max, 100)
        T_points = np.linspace(T_min, T_max, 100)

        # Compute 2D FF array
        FF = np.zeros((len(T_points), len(V_points)))
        P = np.zeros_like(FF)

        for i, T in enumerate(T_points):
            rho_liq = CP.PropsSI("D","T",T,"Q",0,fluid)
            V_liq = m_ox / rho_liq
            P[i,:] = CP.PropsSI("P","T",T,"Q",0,fluid)
            FF[i,:] = V_liq / V_points  # FF for all candidate tank volumes at this T

        # Plot as color map
        fig, ax = plt.subplots(figsize=(8,5))
        c = ax.pcolormesh(V_points, P/1e6, FF, cmap='plasma', shading='auto')  # Pressure in MPa
        # Move colorbar to the very right of the figure
        cbar = fig.colorbar(c, ax=ax, orientation='vertical', fraction=0.015, pad=0.02)
        cbar.set_label("Fill Fraction (FF)")

        ax.set_xlabel("Tank Volume [m^3]")
        ax.set_ylabel("Pressure [MPa]")
        ax.set_title("Nitrous Tank Volume vs Pressure Color Map (Fill Fraction)")

        # Compute pressures for T_min and T_max
        P_min = CP.PropsSI("P","T",T_min,"Q",0,fluid) / 1e6  # MPa
        P_max = CP.PropsSI("P","T",T_max,"Q",0,fluid) / 1e6  # MPa

        # Draw horizontal dashed lines
        ax.hlines([P_min, P_max], V_min, V_max, colors='red', linestyles='dashed')

        # Add labels at the right end of the lines
        ax.text(V_max*1.01, P_min, f'T_min = {T_min-273.15:.1f} °C', va='center', color='red')
        ax.text(V_max*1.01, P_max, f'T_max = {T_max-273.15:.1f} °C', va='center', color='red')

        print("\nUse the graph to estimate the tank volume. Note fill fraction = V_liq/V_tank.\nWe want a tank that isn't too big (long fill time)\nbut we don't want our fill fraction to be too high (overpressurization risk if temp increases)")

        plt.show()


    print(f"\n\n DONE SIZING WIZARD, COPY/PASTE ^ INTO AN INPUT MODULE")





