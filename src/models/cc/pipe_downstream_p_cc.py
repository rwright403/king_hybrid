import pandas as pd
import numpy as np

class pipe_downstream_p_cc:
    def __init__(self, filepath, timestep):
        """
        Parameters
        ----------
        validation_file : str
            Path to CSV with two columns: time, value (P_cc in Pa).
        timestep : float
            Integration timestep [s]
        nozzle : object, optional
            Nozzle model (not required for pure tank coupling)
        P_atm : float, optional
            Ambient pressure [Pa]
        """
        # load validation P_cc data
        df = pd.read_csv(filepath, header=None, skiprows=1, names=["time", "P_cc"])
        self.times = df["time"].values
        self.P_ccs = df["P_cc"].values

        self.t = 0
        self.timestep = timestep

    def inst(self):
        """
        Advance chamber state by one timestep.
        Interpolates P_cc from validation file.
        """
        # interpolate chamber pressure
        self.P_cc = np.interp(self.t, self.times, self.P_ccs)
        self.t += self.timestep

        return {"P_cc": self.P_cc}
