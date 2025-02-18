from scipy.integrate import solve_ivp
import numpy as np
import CoolProp.CoolProp as CP
from thermo.eos import PR
from thermo import Chemical
import matplotlib.pyplot as plt
import traceback

### this is to test:
LOOKUP_TIME = 0

# Global Constants:
R_U = 8.31446 #J/(mol K)
T_REF = 298.15 #K

n2o_g = Chemical('nitrous oxide', T=T_REF)
PC = n2o_g.Pc
TC = n2o_g.Tc
OMEGA = n2o_g.omega

MW = (n2o_g.MW/1000) #n2o.MW in g/mol --> converted to kg/mol
KAPPA = 0.37464 + 1.5422*n2o_g.omega - 0.26992*n2o_g.omega**2
b = 0.07780*(R_U*TC/PC)
g = 9.81 #m/s^2


# Set the conditions and imports
from thermo import ChemicalConstantsPackage, PRMIX, CEOSLiquid, CEOSGas, FlashPureVLS
fluid = 'N2O'
constants, correlations = ChemicalConstantsPackage.from_IDs([fluid])

T = 290
P = 2.5e6
zs = [1]

eos_kwargs = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas)
gas = CEOSGas(PRMIX, HeatCapacityGases=correlations.HeatCapacityGases, 
              eos_kwargs=eos_kwargs)

h = gas.H()/MW
u = gas.U()/MW

print(h, u)
