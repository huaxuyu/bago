Generate search space
---------------------

This function is implemented in the :class:`gpModel` class to generate a search space.

.. function:: genSearchSpace(self, parameters)

    Function to fit the data to Gaussian process regression model.

    :param parameters: dict. A dictionary containing the parameters for finding the top signals.

It generates the following attributes in the :class:`gpModel` object:

    * :attr:`gridX`: a numpy array containing the gradients for the search space.

    * :attr:`gradientPct`: a 1D numpy array containing the percentage of strong mobile phase used by each gradient in the search space.

    * :attr:`scaler`: a :class:`sklearn.preprocessing.StandardScaler` object used to scale the data.

    * :attr:`scaledX`: a numpy array containing the scaled gradients for the search space.

To use this function:

.. code-block::

    # You need a gpModel object (d).

    d.genSearchSpace(parameters)

