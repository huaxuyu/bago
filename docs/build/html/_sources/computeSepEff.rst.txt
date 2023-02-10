Compute Separation Efficiency
-----------------------------

.. function:: sepEfficiency(rtSeq, rtRange)

    Calculate the separation efficiency using a series of retention times.

    :param rtSeq: numpy array. Retention times of the top features.

    :param rtRange: numpy array. Times when gradient begins and ends. In minute.

    :return: float. Separation efficiency.

To use this function:

.. code-block::
    
    sep_eff = sepEfficiency(rtSeq, rtRange)
    print('Separation efficiency: %.2f' % sep_eff)