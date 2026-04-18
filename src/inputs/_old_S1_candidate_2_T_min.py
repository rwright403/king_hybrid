# ------------------------
# Std Engine Input File (Converted Test Case - Liquid Rocket)
# ------------------------
import numpy as np

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
sim_time = 6.7        # [s]
P_atm    = 101325     # [Pa]
T_atm    = 273.15 + 17     # [K]
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
m_ox        =  7.650873122104811   # [kg]
import CoolProp.CoolProp as CP
P_ox_tank   = CP.PropsSI('P', 'T', T_atm, 'Q', 0, 'N2O')          # [Pa]
V_ox_tank      = 0.013   # [m^3]
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
A_inj_ox = 1.5* 3.931764339889601e-05   # [m^2]

# ------------------------
# Fuel tank & pressurant
# ------------------------
m_fuel       = 1.5301746244209624    # [kg]
P_pres_tank  = 7.2e6            # [Pa] fuel tank pressure - used to sol m_pres
V_fuel_tank       = 0.004     # [m^3] fuel tank volume
V_pres_tank  = 0.0             # not provided
diam_out_fuel = 0.0254 * 5.5      # [m]
diam_in_fuel  = 0.0254 * 5.0      # [m]
rho_wall    = 2700          # [kg/m^3]
k_w         = 237           # [W/(m K)]

# ------------------------
# Fuel Injector parameters
# ------------------------
Cd_inj_fuel   = 0.6 #NOTE: GUESS
A_inj_fuel = 0.75*8.353139461101712e-06       # [m^2]

# ------------------------
# Chamber (liquid engine)
# ------------------------
V_cc      = 0.003      # [m^3]

# ------------------------
# Nozzle
# ------------------------
d_throat = 0.0254* 1.35 # [m]
expratio = 6.81741


# ------------------------
# Sizing wizard estimates (kept for reference, not used in runtime)
# ------------------------
min_TW_ratio = 11
Cd_est = 0.6
mass_fraction_estimate = 0.2657741984
characteristic_len = 1.0        #[m]
elevation = 225 # TODO: CHECK ALT MALAHAT RANGE VS TIMMINS
apogee_height = 2*3048      #[m]



# ==========================
# ROCKETPY FLIGHT DESVAR
# ==========================
import datetime
from rocketpy import Environment


# --- LAUNCH CANADA TIMMINS ONT. LAUNCH PAD ---
"""
From Rocketpy Delivarance Mk. IV example flight (Launched at LC 2024)
(no wind, velo components set to 0)
https://docs.rocketpy.org/en/latest/examples/defiance_flight_sim.html
"""
flight_date = datetime.date(2024, 8, 24)
env = Environment(latitude=47.966527, longitude=-81.87413, elevation=295.0)
env.set_date((flight_date.year, flight_date.month, flight_date.day, 0))
env.set_atmospheric_model(type="custom_atmosphere", wind_v=8.3333, wind_u=0.0) # wind_v as 30km/h is the highest wind speed we launch in 

rail_length = 5.6388 #[m] #NOTE: THIS IS NOT THE LC DEFAULT RAIL LENGTH             could possibly adjust this if it solve all our problems, just for learning purpose
inclination = 87     # [degrees, 0 degrees is along ground]                         const.
heading = 90        # [degrees]                                                     const.

# --- Rocket Geometry and Mass Properties ---
fuselage_radius = 0.5*diam_out                  # [m]                               can chamge
fuselage_inner_radius = 0.95*fuselage_radius    #[m]                                match thickness to past year rockets #future work
rkt_csys = "tail_to_nose" # tail = 0 nose = (+) max

upper_launch_lug_pos = 0.8      # [m] csys: "tail_to_nose"                          Ryan will let us know
lower_launch_lug_pos = 0.2      # [m] csys: "tail_to_nose"                          ^
launch_lug_angular_pos = 45     # [degrees]

nose_length = 0.5               # [m]                                               Can change adjust depending on old rockets look at our prior research
nose_kind = "vonKarman"
nose_position = 3.3             # [m]
rho_nose = 2500            # [kg/m^3] assume glass

fins_n = 3                      #NOTE: drag model 2 assumes 4 fin rocket.           keep at 4
fins_span = 0.1016*.9               # [m]                                               base to tip of fin (conctor to the rocket to outside most edge)
fins_root_chord = 0.254*.8           # [m]                                               lostgest distance up    
fins_tip_chord = 0.0762*.8           # [m]                                               shortest distance up
fins_position = 0.3            # [m] csys: "tail_to_nose"                          Change as required "if is hanging off the rocket"
gamma_LE_sweep= 30             # [degrees]
fin_root_thickness=0.0254*(0.25)       # [m]


# fin X-section
Lambda_1= gamma_LE_sweep                 # [degrees]
Lambda_L= Lambda_1              # [degrees]
Lambda_2= 0.0                   # [degrees]
Lambda_T= 0.0                   # [degrees]
zeta_L= 10 #A2 is 30 but this sends C_DWT wayyy too large because Buseman low aoa is violated                 # [degrees]
zeta_T= zeta_L                     # [degrees]
lL_root= 0.0254*0.25            # [m]
lT_root= lL_root                # [m]



rho_fin = 1350 # [kg/m^3] assume cf


ox_tank_pos = 1.1           # [m] csys: "tail_to_nose"
fuel_tank_pos = 2.2
engine_pos = 0              # [m] csys: "tail_to_nose"
cc_cg = 0.15                 # [m] csys: "tail_to_nose"                               
nozzle_pos = 0.0            # [m] csys: "tail_to_nose"

#POINT MASSES OF COMPONENTS:
m_mev = 5.5                 # [kg]                                                                      #man engin center off mass
m_otv = 4.5                 # [kg]
m_reco = 5.5                 # [kg] 
m_ftv= 1.0                 # [kg]

#Densities:
rho_upperfuse = 1350 # [kg/m^3] assume cf
rho_lowerfuse = 1350 # [kg/m^3] assume cf


# mass model = 1 inputs:
#rkt_dry_mass = 37                              # [kg]
#cc_dry_inertia = (1.801, 1.801, 0.0305)        # tuple, [kg m^2]
#cc_dry_mass = 13.832                           # [kg]
#rkt_dry_inertia = (94.14, 94.14, 0.09)         # tuple, [kg m^2]
#rkt_dry_cg = 3.29                              # [m]

# drag model = 1 inputs:
#power_off_drag = "src/inputs/UofT_Deliverance_II/UofT_Style_Defiance_DragCurve.csv"
#power_on_drag = "src/inputs/UofT_Deliverance_II/UofT_Style_Defiance_DragCurve.csv"
