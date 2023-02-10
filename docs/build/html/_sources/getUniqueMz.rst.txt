Get unique m/z
--------------

.. function:: getUniqueMz(d, precsMzTol=0.01)

    Function to get the unique m/z values.

    :param d: a :class:`MSData` object. MS data.

    :param precsMzTol: float. The m/z tolerance to use to get the unique m/z values, in Da.

    :return: List. A list containing the unique m/z values.

To use this function:

.. code-block::

    # You need a MSData object (d).

    mzList = getUniqueMz(d)