Update model
------------

This function is implemented in the :class:`gpModel` class to update the model with new data.

.. function:: updateModel(self, exp, parameters)

    Function to update the model with new data.

    :param exp: dict. A dictionary of LC-MS experiemnts. Key is the name of the experiment and value is a :class:`MSData` object.

    :param parameters: dict. A dictionary of parameters.

To use this function:

.. code-block::

    # You need a gpModel object (model).

    model.updateModel(exp, parameters)