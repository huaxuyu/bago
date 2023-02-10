Generate second gradient
------------------------

.. note::

    When the second gradient is not specified by user, it is generated.

.. function:: computeSecondGradient(parameters, model)

    Calculate the second gradient to run.

    :param parameters: dict. Global parameters.

    :param model: a :class:`gpModel` object. Gaussian process model.

    :return: numpy array. The second gradient.

To use this function:

.. code-block::

    # You need a gpModel object (model)

    computeSecondGradient(parameters, model)