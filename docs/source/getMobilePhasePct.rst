Get mobile phase percentage
---------------------------

.. function:: getMobilePhasePct(gradient, timePoints)

    Function to calculate the percentage of the strong mobile phase used 
    through the gradient, which can be used to estimate the strength of gradient.

    :param gradient: numpy array. Sequence of mobile phase percentages at all time points.

    :param timePoints: numpy array. Sequence of time points.

    :returns: float. The percentage of the strong mobile phase.

To use this function:

.. code-block::

    pctg = getMobilePhasePct(gradient, timePoints)
    print(pctg)