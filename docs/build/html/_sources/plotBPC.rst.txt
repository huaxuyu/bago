Plot BPC
--------

This function is implemented in the :class:`MSData` class to plot the base peak chromatogram (BPC).

.. function:: plotBPC(self, pltName)

    Function to plot BPC.

    :param pltName: str. Name of the plot


To use this function:

.. code-block::

    # You need a MSData object (d)
    
    d.plotBPC('BPC.png')