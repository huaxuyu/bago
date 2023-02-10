MS data (object)
----------------

This is an object for storing, manipulating and visualizing LC-MS data.

The object contains the following attributes:

    * :attr:`rawData`: A :class:`MSExperiment` object (supported by :mod:`pyopenms`).

    * :attr:`ms1Data`: A :class:`ms1Spectrum` object (supported by :mod:`bago`).

    * :attr:`ms2Data`: A :class:`ms2Spectrum` object (supported by :mod:`bago`).

    * :attr:`sNum`: Number of top signals.

    * :attr:`topSignals`: Top signals.

    * :attr:`intensityBPC`: Intensities of the base peak chromatogram.

    * :attr:`rtBPC`: Retention times of the base peak chromatogram.

    * :attr:`sepEff`: Separation efficiency.

Parameters are needed to process the MS data. A templete of the parameters can be accessed by

.. code-block:: python

    from bago import bagoMain

    print(bagoMain.parametersTemplate)