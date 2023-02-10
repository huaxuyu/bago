Gaussian process regression (object)
------------------------------------

This is an object for building Gaussian process regression model and using
it for Bayesian optimization.

The object contains the following attributes:

    * :attr:`kernel`: a :class:`sklearn.gaussian_process.kernels` object.

    * :attr:`gpr`: a :class:`sklearn.gaussian_process.GaussianProcessRegressor` object.

    * :attr:`trainX`: a numpy array containing the training points (input).

    * :attr:`trainy`: a numpy array containing the training values (output).

    * :attr:`scaledX`: a numpy array containing the scaled training points (input).

    * :attr:`gradientPct`: a numpy array containing percentages of the strong mobile phase for gradients in the search space.

    * :attr:`gridX`: a numpy array containing the grid search space.

    * :attr:`scaler`: a :class:`sklearn.preprocessing.StandardScaler` object.