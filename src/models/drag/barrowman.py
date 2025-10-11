import numpy as np
from src.models.drag._base import BaseDragModel



class barrowman(BaseDragModel):
    def __init__(self, geometry=None, method="linear", M_min=0.0, M_max=2.5, n_points=50):
        super().__init__(method=method, M_min=M_min, M_max=M_max, n_points=n_points)
        self.geometry = geometry or {}

    def compute_cd(self, mach_array, x): #this becomes a 2D integration problem

        #NOTE: Vectorized

        # --- Placeholder "fake" Cd curve (replace later with real Barrowman model)
        # Simulates a typical transonic bump
        Cd = 0.3 + 0.01 * np.exp(-((mach_array - 1.0) ** 2) / 0.05)
        Cd += 0.02 * np.sin(2 * np.pi * mach_array / 3.0)

        return Cd