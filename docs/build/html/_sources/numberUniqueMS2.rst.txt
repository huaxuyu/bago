Get unique MS2 spectra
----------------------

.. function:: getUniqueMS2(d, rtTol=1.0, precsMzTol=0.01, dpTol=0.95)

    Function to get the unique MS2 spectra from a :class:`MSData` object.

    :param d: a :class:`MSData` object. MS data.

    :param rtTol: float. Retention time tolerance in minutes.

    :param precsMzTol: float. Precursor m/z tolerance in Da.

    :param dpTol: float. Dot product tolerance.

    :return: List. A list containing the grouped MS2 spectra.

To use this function:

.. code-block::

    # You need a MSData object (d).

    ms2List = getUniqueMS2(d)
