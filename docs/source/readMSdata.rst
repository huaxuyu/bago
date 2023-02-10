Read MS data
------------

This function is implemented in the :class:`MSData` class to read the raw MS data.

mzML and mzXML files are supported. This function is achieved by using the :mod:`pyopenms` module.

.. function:: readRawData(self, fileName)

    Function to read the raw LCMS data to Function to read the raw LCMS data to :class:`MSExperiment` object (supported by :mod:`pyopenms`).

    :param fileName: str. The name of the file to be read (end with .mzML or .mzXML).

The loaded data is stored in the :class:`MSData` object as :class:`MSExperiment` object.
Two more functions are implemented to get the MS1 and MS2 spectra from the :class:`MSExperiment` object.

.. function:: extractMS1(self, rtRange=None)

    Function to extract all MS1 scans and convert them to :class:`ms1Spectrum` objects (supported by bago).

    :param rtRange: list. The retention time range to be extracted (in minute). If None, all MS1 scans will be extracted.

It generate an attribute :attr:`ms1Data` in the :class:`MSData` object to store the extracted MS1 spectra.

.. function:: extractMS2(self, rtRange=None)

    Function to extract all MS2 scans and convert them to :class:`ms2Spectrum` objects (supported by bago).

    :param rtRange: list. The retention time range to be extracted (in minute). If None, all MS2 scans will be extracted.

It generate an attribute :attr:`ms1Data` in the :class:`MSData` object to store the extracted MS1 spectra.

To load the raw data and extract the MS1 and MS2 spectra, you can use the following code:

.. code-block::

    from bago import rawDataHelper

    fileName = 'data.mzML'

    msData = rawDataHelper.MSData()
    msData.readRawData(fileName)
    msData.extractMS1()
    msData.extractMS2()