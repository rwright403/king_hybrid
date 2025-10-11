import numpy as np
from scipy.interpolate import CubicSpline

class BaseDragModel:
    """Base class for all drag models."""

    def __init__(self, method="linear", M_min=0.0, M_max=2.5, n_points=50):
        self.method = method
        self.M_min = M_min
        self.M_max = M_max
        self.n_points = n_points


    def make_drag_model(self, drag_kwargs):
        """Return a callable Cd(M) function suitable for RocketPy."""
        mach_array = np.linspace(self.M_min, self.M_max, self.n_points)
        cd_array = self.compute_cd(mach_array, **drag_kwargs)

        if self.method == "spline":
            interp = CubicSpline(mach_array, cd_array, extrapolate=True)
            def Cd(M):
                return float(interp(M))
        else:
            def Cd(M):
                return float(np.interp(M, mach_array, cd_array))

        Cd.__name__ = self.__class__.__name__
        return Cd

