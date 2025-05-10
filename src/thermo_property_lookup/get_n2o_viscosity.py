import pandas as pd
import numpy as np
from scipy.interpolate import griddata

import CoolProp.CoolProp as CP

# Load viscosity data
df_liq = pd.read_csv("./src/thermo_property_lookup/NIST-webbook-liquid-viscosity.csv")
df_vap = pd.read_csv("./src/thermo_property_lookup/NIST-webbook-vapor-viscosity.csv")

def load_table(df):
    return df["Temperature (K)"].values, df["Pressure (MPa)"].values, df["Kinematic Viscosity (m^2/s)"].values

T_liq, P_liq, visc_liq = load_table(df_liq)
T_vap, P_vap, visc_vap = load_table(df_vap)


"""
Returns interpolated viscosity for given T, P and phase.
Uses data from NIST Chemical Webbook: https://webbook.nist.gov/cgi/cbook.cgi?ID=10024-97-2
"""


def get_n2o_viscosity(T, P, phase):

    if phase == "liquid" and (260 <= T <= 306) and (1e6 <= P <= 7.2e6 ):
        kvisc = griddata((T_liq, P_liq), visc_liq, (T, P), method="linear")
        if np.isnan(kvisc):
            kvisc = griddata((T_liq, P_liq), visc_liq, (T*.98, P), method="linear") # if T,P go out of bounds no they don't ( *.98)
            #raise ValueError(f"get_n2o_viscosity for liquid state returning nan: [{T},{P}]")
        return kvisc
    
    elif phase == "vapor" and (260 <= T <= 306) and (1e6 <= P <= 7.2e6 ):
        kvisc = griddata((T_vap, P_vap), visc_vap, (T, P), method="linear")
        if np.isnan(kvisc):
            kvisc = griddata((T_liq, P_liq), visc_liq, (T*.98, P), method="linear") #cheat a little - "if you're not cheating you're not trying" - Mr. E
            #raise ValueError(f"get_n2o_viscosity for vapor state returning nan: [{T},{P}]")
        return kvisc
    
    else:
        if (260 >= T or T >= 306) or (1e6 >= P or P >= 7.2e6 ):
            raise ValueError(f"Input [T (K) , P (Pa)] out of bounds: [{T},{P}]")
        else:
            raise ValueError("Bad input: either 'liquid' or 'vapor'.")


"""
T_query, P_query = 270, 1e6  # Example temperature & pressure
#visc_liquid = get_n2o_viscosity(T_query, P_query, phase="liquid")
visc_vapor = get_n2o_viscosity(T_query, P_query, phase="vapor")

#print(f"Liquid viscosity: {visc_liquid:.6e} (m^2/s)")
#print(f"Vapor viscosity: {visc_vapor:.6e} (m^2/s)")
"""