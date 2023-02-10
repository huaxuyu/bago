Gaussian process regression fitting
-----------------------------------

The Gaussian process regression fitting is implemented in the
:class:`gpModel` class. You have to generate a search space before
fitting the model.

.. function:: fit(self, X, y)

    Function to fit the data to Gaussian process regression model.

    :param X: numpy array. Training data, input.

    :param y: numpy array. Training data, output.

To use this function:

.. code-block::

    # You need a gpModel object (d) with a search space defined.

    d.fit(X, y)

