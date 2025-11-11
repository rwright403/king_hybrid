# ------------------------
# Std Engine Input File (Converted Test Case - Liquid Rocket)
# ------------------------

# Models
ox_tank_model   = 2   
fuel_tank_model = 1  
ox_inj_model    = 3
cc_model        = 2
nozzle_model    = 1 
drag_model      = 2
mass_model      = 2

# Global settings
thrust_curve_graphs = True
mode = "full_stack"
save_path = None

# Environment
timestep = 0.003       # [s]
sim_time = 8.0        # [s]
P_atm    = 101325     # [Pa]
T_atm    = 293.15     # [K]
rho_atm  = 1.225       # [kg/m^3]

# ------------------------
# Propellant properties
# ------------------------
oxidizer_name   = "N2O"
fuel_name       = "Ethanol"
fuel_properties = None     # not needed for LRE
fuel_str        = "Ethanol"
pres_str        = "N2"     # pressurant gas

# ------------------------
# Oxidizer Tank
# ------------------------
m_ox        = 10.989419425908563   # [kg]
P_ox_tank   = 5.0e6           # [Pa]
V_ox_tank      = 0.0185   # [m^3]
diam_out    = 0.0254*5.5      # [m]
diam_in     = 0.0254*5.0      # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]
volume_err_tol = 0.001
P_dot_err_tol  = None #NOTE: obsolte

# ------------------------
# Ox Injector parameters
# ------------------------
Cd_inj_ox   = 0.6
A_inj_ox = 5.6440e-05   # [m^2]

# ------------------------
# Fuel tank & pressurant
# ------------------------
m_fuel       = 2.747354856477141    # [kg]
m_pres       = 3*0.03428             # [kg] pressurant mass (estimated)
P_pres_tank  = 5.0e6            # [Pa] fuel tank pressure
V_fuel_tank       = 0.006 #0.004     # [m^3] fuel tank volume
V_pres_tank  = 0.0             # not provided
diam_out_fuel = 0.0254 * 5.5      # [m]
diam_in_fuel  = 0.0254 * 5.0      # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]

# ------------------------
# Fuel Injector parameters
# ------------------------
Cd_inj_fuel   = 0.6 #NOTE: GUESS
A_inj_fuel = 1.4977e-05       # [m^2]

# ------------------------
# Chamber (liquid engine)
# ------------------------
V_cc      = 0.003      # [m^3]

# ------------------------
# Nozzle
# ------------------------
d_throat = 0.0254*1.6635 # [m]
expratio = 6.81741


# ------------------------
# Sizing wizard estimates (kept for reference, not used in runtime)
# ------------------------
min_TW_ratio = 11
Cd_est = 0.6
mass_fraction_estimate = 0.2657741984
characteristic_len = 1.0        #[m]
elevation = 225 # TODO: CHECK ALT MALAHAT RANGE VS TIMMINS
apogee_height = 3048*3      #[m]



# ==========================
# ROCKETPY FLIGHT DESVAR
# ==========================
from rocketpy import Environment

# --- LAUNCH CANADA TIMMINS ONT. LAUNCH PAD ---
"""
From Rocketpy Delivarance Mk. IV example flight (Launched at LC 2024)
(no wind, velo components set to 0)
https://docs.rocketpy.org/en/latest/examples/defiance_flight_sim.html
"""

env = Environment(latitude=47.966527, longitude=-81.87413, elevation=1383.4)

rail_length = 5.334 #[m]
inclination = 85
heading = 90

# --- Rocket Geometry and Mass Properties ---
fuselage_radius = 0.07
fuselage_inner_radius = 0.067
rkt_dry_mass = 37.211

upper_launch_lug_pos = 1.25
lower_launch_lug_pos = 0.5
launch_lug_angular_pos = 45


nose_length = 0.563
nose_kind = "vonKarman"
nose_position = 4.947

fins_n = 3
fins_span = 0.115
fins_root_chord = 0.4
fins_tip_chord = 0.2
fins_position = 0.175
gamma_LE_sweep= 15             # [degrees]
fin_root_thickness=1e-3        # [m]

### NOTE: is this cg or location actually??
ox_tank_cg = 2.2
engine_cg = 0
fuel_tank_cg = 1.0

# --- Chamber / Motor Properties ---
cc_cg = 0.25
cc_dry_inertia = (1.801, 1.801, 0.0305)
cc_dry_mass = 13.832
nozzle_pos = 0.0

"""
############
Preliminary Design Summary:
############


Performance:
------------
Specific Impulse: 270.15145811600627 (s)
Characteristic Velocity: 1578.597237970397 (m/s)
Burn Time: 4.5 (s)
Rocket Dry Mass: 51.6858836 (kg)
Total Impulse: 29792.12 (N s)

Combustion Chamber:
------------
P_cc: 3000000.0 (Pa)
O/F Ratio: 4.0
Flame Temperature: 3065.766452364495
Ratio of Specific Heats: 1.1680264089309822
Reactant Gas Const. 336.7284307249569 (J/(kg K))
Average Thrust: 6499.171461737702 (N)
Mass Flow Rate: 2.555611631851842 (kg/s)
Exit Velocity: 2543.098247298353 (m/s)
Expansion Ratio: 6.817409006311842
Throat Diameter 1.663535923934223 (in)

Feed System:
------------
Estimated Mass Fraction: 0.2657741984
Fuel Mass: 2.747354856477141 (kg)
Oxidizer Mass: 10.989419425908563 (kg)
Oxidizer Tank Pressure: 5000000.0 (Pa)
Fuel Tank Pressure: 5000000.0 (Pa)
Preliminary Injector Solved by Assuming Fuel and Oxidizer Orifices Follow SPI Model and a Cd guess of 0.6
Oxidizer Injector Discharge Area: 5.644023527701955e-05 (m^2), Fuel Injector Discharge Area 1.4977972405728426e-05 (m^2)

fuel tank summary: m_pressurant = 0.03427696559234706, V_tank = 0.0039968207632066505
"""