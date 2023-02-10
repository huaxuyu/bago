Output gradient
---------------

.. function:: outputConfig(name, timePoints, gradient)

    Function to output the gradient to a csv file.

    :param name: str. Name of the file to be output.

    :param timePoints: numpy array. Time points for the gradient.

    :param gradient: numpy array. Gradient to be output.

To use this function:

.. code-block::

    from bago import rawDataHelper
    rawDataHelper.outputConfig('gradient.csv', timePoints, gradient)