Get BPC
-------

This function is implemented in the :class:`MSData` class to obtain the base peak chromatogram (BPC).

.. function:: getBPCData(self)

    Function to prepare data for BPC.

It generates the following attributes in the :class:`MSData` object:

    * :attr:`rtBPC`: retention time points of the BPC
    * :attr:`intensityBPC`: intensity values of the BPC

To use this function:

.. code-block::

    # You need a MSData object (d)

    d.getBPCData()
