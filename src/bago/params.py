# Author: Huaxu Yu

import numpy as np
import pandas as pd
import os
import pickle

from .raw_data_utils import read_method_from_csv


class Parameters:

    def __init__(self):

        self.time_arr: np.ndarray = None
        self.fix_time2_pct: bool = False
        self.pct_range: tuple = None
        self.pct_step: int = None
        self.total_mobile_phase_ratio_range: tuple = None   # range of the total mobile phase ratio
        self.project_dir: str = None

        # default values recommended
        self.gaussian_sim_tol: float = 0.8                  # tolerance for the Gaussian similarity
        self.noise_score_tol: float = 0.3                   # tolerance for the noise score
        self.acq_func: str = 'ei'                           # acquisition function

        # no need to set these parameters
        self.experiments: list = []
        self.gradients: dict = {}
        self.valid_pcts: np.ndarray = None
        self.rt_range: tuple = None
        self.data_dir: str = None
        self.method_dir: str = None


    def prepare_params(self):
        """
        Prepare the parameters.
        """

        self.data_dir = os.path.join(self.project_dir, "data")
        self.method_dir = os.path.join(self.project_dir, "method")

        # read the time array from the first method file
        path = os.path.join(self.method_dir, "1.csv")
        if not os.path.exists(path):
            raise ValueError("The first method file does not exist. Please check the method directory.")
        else:
            self.time_arr, _ = read_method_from_csv(path)
        
        self.rt_range = (self.time_arr[0], self.time_arr[-1])

        p1 = int(self.pct_range[0])
        p2 = int(self.pct_range[1])
        self.valid_pcts = np.arange(p1, p2, self.pct_step, dtype=int)
        self.valid_pcts = np.append(self.valid_pcts, self.pct_range[1])


    def save_params(self, dir=None):
        """
        Save the parameters to a file.
        """

        # Save the parameters
        if dir is not None:
            pickle.dump(self, open(os.path.join(dir, "params.pkl"), "wb"))
        elif self.project_dir is not None:
            pickle.dump(self, open(os.path.join(self.project_dir, "params.pkl"), "wb"))
        else:
            raise ValueError("Project directory is not set. Please set the project directory before saving the parameters.")