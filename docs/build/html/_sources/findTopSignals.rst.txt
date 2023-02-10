Find top signals
------------------

A function in the :class:`MSData` class to find the top signals in the MS1. 

.. function:: findTopSignals(self, parameters)

    Function to find the top signals in the raw LCMS data.

    :param parameters: dict. A dictionary containing the parameters for finding the top signals.

To use this function:

.. code-block::

    # You need a MSData object (d).

    d.findTopSignals(parameters)