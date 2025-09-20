import cantera as ct
import numpy as np
from scipy.interpolate import RegularGridInterpolator


import cantera as ct
from rocketcea.cea_obj_w_units import CEA_Obj

def build_cantera_from_cea(C, oxName="N2O", fuelName="C2H5OH",Pc=3e6, OF=6.0, eps=40.0,nasa_file="nasa_gas.yaml"):
    #TODO:

    """
    Use RocketCEA to determine relevant species, then
    build a Cantera Solution with those species.

    Parameters
    ----------
    oxName : str
        RocketCEA oxidizer name
    fuelName : str
        RocketCEA fuel name
    Pc : float
        Chamber pressure [Pa]
    OF : float
        Oxidizer/fuel mass ratio
    eps : float
        Expansion ratio
    nasa_file : str
        Cantera NASA species database yaml

    Returns
    -------
    gas : ct.Solution
        Cantera gas object with only the relevant species
    rocketcea_species : list[str]
        List of species names (RocketCEA → Cantera mapping)
    """

    # get equilibrium species set from RocketCEA
    molWtD, massFracD = C.get_SpeciesMassFractions(Pc=Pc, MR=OF, eps=eps)

    # load full NASA species database
    full_species = {S.name: S for S in ct.Species.list_from_file(nasa_file)}

    # map RocketCEA → Cantera names
    def map_species_name(name):
        return name.replace("*", "")  # "*O" -> "O", "*N2" -> "N2", etc.

    rocketcea_species = [map_species_name(sp) for sp in massFracD.keys()]

    # ensure reactants are included
    rocketcea_species += [fuelName, oxName]

    # filter to those available in Cantera’s DB
    species = [full_species[s] for s in set(rocketcea_species) if s in full_species]

    # build minimal Cantera phase
    gas = ct.Solution(thermo="ideal-gas", kinetics="none", species=species)

    return gas, rocketcea_species


import cantera as ct
import numpy as np
from scipy.interpolate import RegularGridInterpolator

class CombustionChamberLookup:
    def __init__(self,C, fuel="C2H5OH", oxidizer="N2O",
                 rho_range=(0.1, 40.0),     # kg/m^3
                 u_range=(-2e6, 5e6),       # J/kg (specific internal energy)
                 OF_range=(0.5, 10.0),      # O/F ratios
                 nRho=40, nU=100, nOF=20):

        # Build a Cantera gas object from RocketCEA species list
        self.gas, sp_list = build_cantera_from_cea(C, oxidizer, fuel, Pc=3e6, OF=6.0)

        # grid axes
        self.rho_grid = np.linspace(*rho_range, nRho)
        self.u_grid   = np.linspace(*u_range, nU)
        self.OF_grid  = np.linspace(*OF_range, nOF)

        # 3D tables (rho × OF × u)
        shape = (nRho, nOF, nU)
        self.T_table     = np.zeros(shape)
        self.P_table     = np.zeros(shape)
        self.h_table     = np.zeros(shape)
        self.cp_table    = np.zeros(shape)
        self.gamma_table = np.zeros(shape)
        self.MW_table    = np.zeros(shape)

        MW_fuel = self.gas.molecular_weights[self.gas.species_index(fuel)]
        MW_ox   = self.gas.molecular_weights[self.gas.species_index(oxidizer)]
        R_univ  = ct.gas_constant

        # Sweep grid
        for i, rho in enumerate(self.rho_grid):
            for j, OF in enumerate(self.OF_grid):
                # Convert OF → mole fractions
                n_fuel = 1.0
                m_fuel = n_fuel * MW_fuel
                m_ox   = OF * m_fuel
                n_ox   = m_ox / MW_ox
                X = {fuel: n_fuel, oxidizer: n_ox}

                for k, u_target in enumerate(self.u_grid):
                    # Start with a temperature guess and ideal gas P guess
                    T_guess = 2500.0
                    P_guess = rho * R_univ / self.gas.mean_molecular_weight * T_guess

                    try:
                        # Set enthalpy + P + composition
                        h_target = u_target + R_univ * T_guess / self.gas.mean_molecular_weight
                        self.gas.HPX = h_target, P_guess, X
                        self.gas.equilibrate("HP")

                        # Store properties
                        self.T_table[i,j,k]     = self.gas.T
                        self.P_table[i,j,k]     = self.gas.P
                        self.h_table[i,j,k]     = self.gas.enthalpy_mass
                        self.cp_table[i,j,k]    = self.gas.cp_mass
                        self.gamma_table[i,j,k] = self.gas.cp_mass / self.gas.cv_mass
                        self.MW_table[i,j,k]    = self.gas.mean_molecular_weight
                    except Exception:
                        # If equilibrium fails, store NaNs
                        self.T_table[i,j,k]     = np.nan
                        self.P_table[i,j,k]     = np.nan
                        self.h_table[i,j,k]     = np.nan
                        self.cp_table[i,j,k]    = np.nan
                        self.gamma_table[i,j,k] = np.nan
                        self.MW_table[i,j,k]    = np.nan

        # Build interpolators on (rho, OF, u)
        self.grid_axes = (self.rho_grid, self.OF_grid, self.u_grid)

        self.T_interp     = RegularGridInterpolator(self.grid_axes, self.T_table, bounds_error=False, fill_value=None)
        self.P_interp     = RegularGridInterpolator(self.grid_axes, self.P_table, bounds_error=False, fill_value=None)
        self.h_interp     = RegularGridInterpolator(self.grid_axes, self.h_table, bounds_error=False, fill_value=None)
        self.cp_interp    = RegularGridInterpolator(self.grid_axes, self.cp_table, bounds_error=False, fill_value=None)
        self.gamma_interp = RegularGridInterpolator(self.grid_axes, self.gamma_table, bounds_error=False, fill_value=None)
        self.MW_interp    = RegularGridInterpolator(self.grid_axes, self.MW_table, bounds_error=False, fill_value=None)

    def query(self, u, rho, OF):
        """
        Query thermo properties at given (u [J/kg], rho [kg/m3], OF).
        Returns T, P, h, cp, gamma, MW.
        """
        pt = (rho, OF, u)
        props = {
            "T":     float(self.T_interp(pt)),
            "P":     float(self.P_interp(pt)),
            "h":     float(self.h_interp(pt)),
            "cp":    float(self.cp_interp(pt)),
            "gamma": float(self.gamma_interp(pt)),
            "MW":    float(self.MW_interp(pt)),
        }
        #               "T"                         "P"                     "h"                         "cp"                      "gamma"                   "MW"
        return float(self.T_interp(pt)), float(self.P_interp(pt)), float(self.h_interp(pt)), float(self.cp_interp(pt)), float(self.gamma_interp(pt)), float(self.MW_interp(pt))
        """NOTE: returns props in order of: T, P, h, cp, gamma, MW"""


import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

def check_monotonicity_hp(fuel="C2H5OH", oxidizer="N2O",
                          P=3e6, OF=4.0, h_min=1.2e6, h_max=1.8e6, nH=200):
    """
    Check whether u(h) is monotonic at fixed pressure and O/F ratio
    using HP (constant enthalpy + pressure) equilibrium.
    """
    gas, sp_list = build_cantera_from_cea(oxidizer, fuel, Pc=P, OF=OF)

    # Build mixture composition from O/F
    MW_fuel = gas.molecular_weights[gas.species_index(fuel)]
    MW_ox   = gas.molecular_weights[gas.species_index(oxidizer)]

    n_fuel = 1.0
    m_fuel = n_fuel * MW_fuel
    m_ox   = OF * m_fuel
    n_ox   = m_ox / MW_ox
    X = {fuel: n_fuel, oxidizer: n_ox}

    # Sweep enthalpy values
    h_range = np.linspace(h_min, h_max, nH)
    u_vals, T_vals = [], []

    for h in h_range:
        try:
            gas.HPX = h, P, X
            gas.equilibrate("HP")
            u_vals.append(gas.int_energy_mass)
            T_vals.append(gas.T)
        except Exception:
            u_vals.append(np.nan)
            T_vals.append(np.nan)

    u_vals = np.array(u_vals)
    T_vals = np.array(T_vals)

    # Check monotonicity
    du = np.diff(u_vals)
    monotonic = np.all(du > 0) or np.all(du < 0)

    # Plot
    fig, ax1 = plt.subplots()
    ax1.plot(h_range, u_vals, label="u(h)")
    ax1.set_xlabel("h [J/kg]")
    ax1.set_ylabel("u [J/kg]")
    ax1.set_title(f"O/F={OF}, P={P/1e5:.1f} bar\nMonotonic: {monotonic}")
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(h_range, T_vals, "r--", label="T(h)")
    ax2.set_ylabel("T [K]")
    plt.show()

    return monotonic, h_range, u_vals, T_vals



#heck_monotonicity(fuel="C2H5OH", oxidizer="N2O",P=3e6, OF=4.0, T_min=500, T_max=4000, nT=200)
check_monotonicity_hp(P=6e6, OF=10) #NOTE: yes, swept a range from (P=.5e6, OF=0.5) to (P=6e6, OF=10), all monotonic! 