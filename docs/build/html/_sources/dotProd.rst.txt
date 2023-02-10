Spectral similarity
--------------------

.. function:: dotProd(ms2A, ms2B, mzTol=0.02)

    Calculate the MS2 spectral similarity by dot product.

    :param ms2A: a :class:`ms2Spectrum` object. The first MS2. 

    :param ms2B: a :class:`ms2Spectrum` object. The second MS2.

    :param mzTol: float. Tolerance of the m/z difference in Da. Default is 0.02.
        
    :return: float. The dot product similarity score between the two MS2 spectra. Returns 2.0 if either MS2 does not exist.

To use this function:

.. code-block::
    
    dp = dotProd(ms2A, ms2B)
    print(dp)