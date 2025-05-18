# imports
import numpy as np
from scipy.stats import norm
from math import pi
from itertools import combinations_with_replacement

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
from sklearn.preprocessing import StandardScaler

from .raw_data_utils import cal_mobile_phase_ratio


class gpModel:
    """
    Class represents the Gaussian process regression model.
    """

    def __init__(self):
        """
        Function to initiate the Gaussian process model by 
        defining the kernel function.
        ----------------------------------------------------------
        """

        # Initialize the kernel function.
        self.kernel = (
                2.0 * Matern(nu=1.5)
        )

        self.gpr = GaussianProcessRegressor(kernel=self.kernel)
        
        self.gradient_arrs = None   # search space (mobile phase pct only at tunable time points)
        self.scaler = None          # StandardScaler object for scaling the search space
        self.scaledX = None         # search space after scaling
        self.total_ratios = None    # total mobile phase ratios
        self.start_pos = None       # start position of the tunable gradient percentage

        self.trainX = None          # training data, X
        self.trainy = None          # training data, y   


    def transform_and_fit(self):
        """
        Function to fit the data to Gaussian process regression model.
        """

        self.gpr.fit(self.scaler.transform(self.trainX[:, self.start_pos:-1]), self.trainy)


    def generate_search_space(self, parameters):
        """
        Generate a search space with all reasonable gradient programs.

        Parameters
        ----------------------------------------------------------
        parameters: dict
            Dictionary of parameters.
        """

        if parameters.fix_time2_pct:
            num = len(parameters.time_arr) - 3
        else:
            num = len(parameters.time_arr) - 2

        # Generate the search space
        gradients = []

        for g in combinations_with_replacement(parameters.valid_pcts, num):
            gradients.append(g)
        
        gradients = np.array(gradients)
        
        a = np.ones((len(gradients), 1), dtype=int) * parameters.pct_range[0]
        c = np.ones((len(gradients), 1), dtype=int) * parameters.pct_range[1]
        
        if parameters.fix_time2_pct:
            b = np.ones((len(gradients), 1), dtype=int) * parameters.pct_range[0]
            self.gradient_arrs = np.concatenate((a, b, gradients, c), axis=1)
            self.start_pos = 2
        else:
            self.gradient_arrs = np.concatenate((a, gradients, c), axis=1)
            self.start_pos = 1
        
        ratios = cal_mobile_phase_ratio(parameters.time_arr, self.gradient_arrs)
        v = (ratios > parameters.total_mobile_phase_ratio_range[0]) & (ratios < parameters.total_mobile_phase_ratio_range[1])

        self.gradient_arrs = self.gradient_arrs[v]

        self.total_ratios = ratios[v]

        # Scale the search space
        self.scaler = StandardScaler()
        self.scaledX = self.scaler.fit_transform(self.gradient_arrs[:, self.start_pos:-1])


    def compute_next_gradient(self, acq="ei"):
        """
        Function to calculate the next gradient to run using an
        acquisition function.

        Parameters
        ----------------------------------------------------------
        acq: str
            Acquisition functions in lower case.
            "ei": expected improvement
            "explore": pure exploration
            "exploit": pure exploitation
            "eps": epsilon-greedy algorithm
            "pi": probability of improvement
            "rand": random a gradient
            "ucb": upper confidence bound
            ""

        Returns
        ----------------------------------------------------------
        numpy array:
            next gradient to run
        """

        best_y = np.max(self.trainy)
        idx = None

        if acq.lower() == "ei":
            idx = expected_improvement(model=self, best_y=best_y)
        elif acq.lower() == "pi":
            idx = probability_of_improvement(model=self, best_y=best_y)
        elif acq.lower() == "eps":
            idx = epsilon_greedy(model=self)
        elif acq.lower() == "exploit":
            idx = exploit(model=self)
        elif acq.lower() == "explore":
            idx = explore(model=self)
        elif acq.lower() == "rand":
            idx = random_gradient(model=self)

        return self.gradient_arrs[idx]


    def update_model(self, parameters):
        """
        Function to update the model by training the model with
        the new data.

        Parameters
        ----------------------------------------------------------
        parameters: Parameters object
            Get the new data from paarameters.experiments. 
        """

        self.trainX = np.array([e.pct_arr for e in parameters.experiments])
        self.trainy = np.array([e.good_feat_num for e in parameters.experiments])
        self.transform_and_fit()


"""
============================================================================================
Acquisition functions
============================================================================================
"""

def expected_improvement(model: gpModel, best_y: float, jitter: float = 0.01):
    """
    Compute expected improvement.

    Parameters
    ----------
    model: gpModel
        Gaussian process regression model.
    best_y: float
        Maximal y observed.
    jitter : float
        Parameter which controls the degree of exploration.

    Returns
    -------
    idx: int
        Index of the selected gradient.
    """

    mean, stdev = model.gpr.predict(model.scaledX, return_std=True)
    stdev = stdev + 1e-6

    # EI parameter values
    z = (mean - best_y - jitter) / stdev
    imp = mean - best_y - jitter
    ei = imp * norm.cdf(z) + stdev * norm.pdf(z)

    return np.argmax(ei)


def probability_of_improvement(model: gpModel, best_y: float, jitter: float = 0.01):
    """
    Compute probability of improvement.

    Parameters
    ----------
    model: gpModel
        Gaussian process regression model.
    best_y: float
        Maximal y observed.
    jitter : float
        Parameter which controls the degree of exploration.

    Returns
    -------
    idx: int
        Index of the selected gradient.
    """

    # Mean and standard deviation
    mean, stdev = model.gpr.predict(model.scaledX, return_std=True)

    # PI parameter values
    z = (mean - best_y - jitter) / stdev
    cdf = norm.cdf(z)

    return np.argmax(cdf)


def epsilon_greedy(model: gpModel, eps: float = 0.1, state: int = 419):
    """
    Epsilon greedy algorithm.

    Parameters
    ----------
    model: gpModel
        Gaussian process regression model.
    eps: float
        Probability of exploration.

    Returns
    -------
    idx: int
        Index of the selected gradient.
    """

    np.random.seed(state)
    
    p = np.random.random()
    if p < eps:
        return random_gradient(model)
    else:
        return exploit(model)


def exploit(model: gpModel):
    """
    Pure exploitation.
    
    Parameters
    ----------
    model: gpModel
        Gaussian process regression model.

    Returns
    -------
    idx: int
        Index of the selected gradient.
    """

    mean = model.gpr.predict(model.scaledX)
    idx = np.argmax(mean)
    # make sure the max is not the one that has been already tested
    a = model.gradient_arrs[idx]
    b = model.trainX
    while np.sum(np.all(b == a, axis=1) > 0) > 0:
        idx = idx + 1
        a = model.gradient_arrs[idx]
    
    return idx


def explore(model: gpModel):
    """
    Pure exploration.

    Parameters
    ----------
    model: gpModel
        Gaussian process regression model.

    Returns
    -------
    idx: int
        Index of the selected gradient.
    """

    _, stdev = model.gpr.predict(model.scaledX, return_std=True)

    return np.argmax(stdev)


def random_gradient(model: gpModel, state: int = 419):
    """
    Randomly select gradient.

    Parameters
    ----------
    model: gpModel
        Gaussian process regression model.
    state: int
        Random state.

    Returns
    -------
    idx: int
        Index of the selected gradient.
    """

    np.random.seed(state)

    return np.random.randint(len(model.gradient_arrs))


def compute_the_second_gradient(parameters, model, ratio_diff=0.1):
    """
    Calculate the second gradient to run.

    Parameters
    ----------
    parameters: Parameters object
        Parameters.
    model: gpModel
        Gaussian process regression model.
    fac: float
        Factor to calculate the second gradient.
    """

    if len(parameters.experiments) < 1:
        raise ValueError("The second gradient can only be calculated when the first gradient is provided.")
    elif len(parameters.experiments) > 1:
        raise ValueError("The second gradient has been provided and threfore no need to generate again.")
    else:
        ratio1 = parameters.experiments[0].total_mobile_phase_ratio
    
    for idx in range(len(model.total_ratios)):
        if abs(ratio1 - model.total_ratios[idx]) > ratio_diff:
            break

    return model.gradient_arrs[idx]