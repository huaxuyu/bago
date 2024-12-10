# imports
import numpy as np
from scipy.stats import norm
from math import pi
from itertools import combinations_with_replacement
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
from sklearn.preprocessing import StandardScaler

from . import rawDataHelper


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

        # Initialize the Gaussian process regression model. 
        # Pass an int for reproducible results across multiple function calls.
        self.gpr = GaussianProcessRegressor(kernel=self.kernel)

        # Search space after scaling.
        self.scaledX = None
        # Mobile phase percentages of the search space.
        self.gradientPct = None
        # Training data.
        self.trainX = None
        self.trainy = None
        # Scaler for the search space.
        self.scaler = None
        # Search space before scaling.
        self.gridX = None


    def fit(self, X, y):
        """
        Function to fit the data to Gaussian process regression model.

        Parameters
        ----------------------------------------------------------
        X: numpy.ndarray
            Training data, X
        y: numpy array
            Training data, y.
        """

        self.trainX = X
        self.trainy = y
        self.gpr.fit(self.scaler.transform(X), y)


    def genSearchSpace(self, parameters):
        """
        Function to generate search space.

        Parameters
        ----------------------------------------------------------
        parameters: dict
            Dictionary of parameters.
        """

        # Calculate the feature number
        fNumber = np.count_nonzero(parameters["isChangable"])

        # Generate the search space
        gradients = []
        temp = np.copy(parameters['grads']['Init_1'])
        gradientPct = []
        for g in combinations_with_replacement(parameters['gradPoints'], fNumber):
            temp[parameters["isChangable"]] = g
            # Calculate the mobile phase percentage
            pct = rawDataHelper.getMobilePhasePct(temp, parameters['timePoints'])
            # Check if the mobile phase percentage is in the range
            if parameters['mpBound'][0] < pct < parameters['mpBound'][1]:
                gradients.append(g)
                gradientPct.append(pct)
        self.gridX = np.array(gradients)
        self.gradientPct = np.array(gradientPct)

        # Scale the search space
        self.scaler = StandardScaler()
        self.scaledX = self.scaler.fit_transform(self.gridX)


    def computeNextGradient(self, acqFunc="ei"):
        """
        Function to calculate the next gradient to run using an
        acquisition function.

        Parameters
        ----------------------------------------------------------
        acqFunc: str
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

        maxObserved = np.amax(self.trainy)
        idx = -1
        if acqFunc.lower() == "ei":
            idx = expectedImprovement(self.gpr, self.scaledX, maxObserved)
        elif acqFunc.lower() == "explore":
            idx = explore(self.gpr, self.scaledX)
        elif acqFunc.lower() == "rand":
            idx = randX(self.gridX, self.trainX)
        elif acqFunc.lower() == "eps":
            idx = epsilonGreedy(self.gpr, self.gridX, self.scaledX, self.trainX)
        elif acqFunc.lower() == "exploit":
            idx = exploit(self.gpr, self.gridX, self.scaledX, self.trainX)
        elif acqFunc.lower() == "pi":
            idx = probabilityOfImprovement(self.gpr, self.scaledX, maxObserved)
        elif acqFunc.lower() == "ucb":
            idx = ucb(self.gpr, self.scaledX, len(self.trainX))
        if idx != -1:
            return self.gridX[idx]
    

    def updateModel(self, exp, parameters):
        """
        Function to update the model by training the model with
        the new data.

        Parameters
        ----------------------------------------------------------
        exp: dict
            Dictionary of the experiment.
        parameters: dict
            Global parameters.
        """

        X = np.array(list(parameters['grads'].values()))
        X = X[:,parameters['isChangable']]
        temp = parameters['grads'].keys()
        y = np.array([exp[n].sepEff for n in temp])
        self.fit(X, y)


def expectedImprovement(model, X, maxObserved, jitter=0.01):
    """
    Compute expected improvement.
    
    EI attempts to balance exploration and exploitation by accounting
    for the amount of improvement over the best observed value.

    Parameters
    ----------------------------------------------------------
    model: sklearn.gaussian_process._gpr.GaussianProcessRegressor
        Trained model.
    X: numpy.ndarray
        Search space.
    maxObserved: float
        Maximal y observed.
    jitter : float
        Parameter which controls the degree of exploration.

    Returns
    ----------------------------------------------------------
    idx: int
        Index of the selected gradient.
    """

    mean, stdev = model.predict(X, return_std=True)
    stdev = stdev + 1e-6

    # EI parameter values
    z = (mean - maxObserved - jitter) / stdev
    imp = mean - maxObserved - jitter
    ei = imp * norm.cdf(z) + stdev * norm.pdf(z)

    return np.argmax(ei)


def probabilityOfImprovement(model, X, maxObserved, jitter=0.01):
    """
    Compute expected improvement.
    
    EI attempts to balance exploration and exploitation by accounting
    for the amount of improvement over the best observed value.

    Parameters
    ----------------------------------------------------------
    model: sklearn.gaussian_process._gpr.GaussianProcessRegressor
        Trained model.
    X: numpy.ndarray
        Search space.
    maxObserved: float
        Maximal y observed.
    jitter : float
        Parameter which controls the degree of exploration.

    Returns
    ----------------------------------------------------------
    idx: int
        Index of the selected gradient.
    """

    # Mean and standard deviation
    mean, stdev = model.predict(X, return_std=True)

    # PI parameter values
    z = (mean - maxObserved - jitter) / stdev
    cdf = norm.cdf(z)

    return np.argmax(cdf)


def explore(model, X):
    """
    Find the gradient with the largest variance

    Parameters
    ----------------------------------------------------------
    model: sklearn.gaussian_process._gpr.GaussianProcessRegressor
        Trained model.
    X: numpy.ndarray
        Search space.

    Returns
    ----------------------------------------------------------
    idx: int
        Index of the selected gradient.
    """

    _, stdev = model.predict(X, return_std=True)

    return np.argmax(stdev)


def randX(gridX, trainX):
    """
    Randomly select gradient.

    Parameters
    ----------------------------------------------------------
    gridX: numpy.ndarray
        gridX.
    trainX: numpy.ndarray
        Training data

    Returns
    ----------------------------------------------------------
    idx: int
        Index of the selected gradient.
    """

    idx = np.random.randint(len(gridX))
    while checkDuplicates(gridX[idx], trainX):
        idx = np.random.randint(len(gridX))
    return idx


def epsilonGreedy(model, gridX, scaledX, trainX, eps=0.05):
    """
    Acquisition function by epsilon-greedy algorithm.
    
    Parameters
    ----------------------------------------------------------
    model: sklearn.gaussian_process._gpr.GaussianProcessRegressor
        Trained model.
    gridX: numpy.ndarray
        gridX.
    scaledX: numpy.ndarray
        scaledX.
    trainX: numpy.ndarray
        Training data
    eps: float
        Probability for random exploration.

    Returns
    ----------------------------------------------------------
    idx: int
        Index of the selected gradient.
    """

    p = np.random.random()
    if p < eps:
        return randX(gridX, trainX)
    else:
        return exploit(model, gridX, scaledX, trainX)


def exploit(model, gridX, scaledX, trainX):
    """
    Acquisition function by epsilon-greedy algorithm.
    
    Parameters
    ----------------------------------------------------------
    model: sklearn.gaussian_process._gpr.GaussianProcessRegressor
        Trained model.
    gridX: numpy.ndarray
        gridX.
    scaledX: numpy.ndarray
        scaledX.
    trainX: numpy.ndarray
        Training data.

    Returns
    ----------------------------------------------------------
    idx: int
        Index of the selected gradient.
    """

    mean = model.predict(scaledX)
    counter = 0
    t = (-mean).argsort()[counter]
    while checkDuplicates(gridX[t], trainX):
        counter += 1
        t = (-mean).argsort()[counter]

    return t


def ucb(model, X, iters):
    """
    Acquisition function by upper confidence bound algorithm.
    
    Parameters
    ----------------------------------------------------------
    model: sklearn.gaussian_process._gpr.GaussianProcessRegressor
        Trained model.
    X: numpy.ndarray
        Search space.
    iters: int
        Iteration number.

    Returns
    ----------------------------------------------------------
    idx: int
        Index of the selected gradient.
    """

    # Mean and standard deviation
    mean, stdev = model.predict(X, return_std=True)

    d = len(X[0])
    coeff = np.sqrt(2*np.log(iters**(d/2+2) * pi**2/3))
    ucb = mean + coeff * stdev

    return np.argmax(ucb)


def checkDuplicates(arr, X):
    """
    A function to check if a numpy array is in a 2D numpy array.
    
    Parameters
    ----------------------------------------------------------
    arr: numpy array
    X: numpy.ndarray

    Returns
    ----------------------------------------------------------
    boolean
        True for in.
    """
    for i in X:
        if np.array_equal(i, arr):
            return True

    return False