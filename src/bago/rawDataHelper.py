from pyopenms import MSExperiment, MzMLFile, MzXMLFile
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt


class MSData:
    """
    A class that models the raw LCMS data in a single file (mzML or mzXML).
    It provides functions to process the raw data.
    """

    def __init__(self):
        """
        Function to initiate MSData.
        ----------------------------------------------------------
        """

        # Raw LCMS data
        self.rawData = None
        # MS1 data
        self.ms1Data = None
        # MS2 data
        self.ms2Data = None

        # Top signals
        self.sNum = None    # Number of top signals
        self.topSignals = None    # Top signals

        self.intensityBPC = None
        self.rtBPC = None
        self.sepEff = None
        

    def readRawData(self, fileName):
        """
        Function to read the raw LCMS data to MSExperiment object 
        (supported by pyopenms package).

        Parameters
        ----------------------------------------------------------
        fileName: str
            File name to read the raw LCMS data, which can be either mzML
            or mzXML (default).
        """

        if os.path.isfile(fileName):
            self.rawData = MSExperiment()
            # get extension from file name
            ext = os.path.splitext(fileName)[1]

            if ext.lower() == ".mzxml":
                MzXMLFile().load(fileName, self.rawData)
            elif ext.lower() == ".mzml":
                MzMLFile().load(fileName, self.rawData)
            else:
                print("Unsupported raw LCMS data format.")
        else:
            print("File does not exist.")


    def clearRawData(self):
        """
        Function to remove the raw LCMS data.
        ----------------------------------------------------------
        """

        self.rawData = None


    def extractMS1(self, rtRange=None):
        """
        Function to extract all MS1 scans and convert them to 
        ms1Spectrum objects.

        Parameters
        ----------------------------------------------------------
        rtRange: list
            A list with two numeric values specifying the range of 
            retention time that the MS1 scans are extracted, in minute.
        """

        # If rtRange is not specified, extract all MS1 scans
        if rtRange is None:
            rtRange = (0, np.inf)
        
        # Extract MS1 scans
        self.ms1Data = []
        for spec in self.rawData:
            # Check if the retention time is within the range
            if rtRange[0] * 60 < spec.getRT() < rtRange[1] * 60:
                # Check if the scan is MS1 and the scan is not empty
                if spec.getMSLevel() == 1 and len(spec.get_peaks()[0]) != 0:
                    # Get m/z, intensity and retention time
                    mz = spec.get_peaks()[0]
                    intensity = spec.get_peaks()[1]
                    rt = spec.getRT() / 60
                    # Create a ms1Spectrum object
                    temp = ms1Spectrum(mz=mz, intensity=intensity, rt=rt)
                    self.ms1Data.append(temp)


    def extractMS2(self, rtRange=None):
        """
        Function to extract all MS2 scans and convert them to 
        ms2Spectrum objects.

        Parameters
        ----------------------------------------------------------
        rtRange: list
            A list with two numeric values specifying the range of
            retention time that the MS2 scans are extracted, in minute.
        """

        # If rtRange is not specified, extract all MS2 scans
        if rtRange is None:
            rtRange = (0, np.inf)

        # Extract MS2 scans
        self.ms2Data = []
        for spec in self.rawData:
            # Check if the retention time is within the range
            if rtRange[0] * 60 < spec.getRT() < rtRange[1] * 60:
                # Check if the scan is MS2 and the scan is not empty
                if spec.getMSLevel() == 2 and len(spec.get_peaks()[0]) != 0:
                    # Get precursor m/z and retention time
                    precsMz = spec.getPrecursors()[0].getMZ()
                    rt = spec.getRT() / 60
                    # Get product m/z and intensity
                    prodMz = spec.get_peaks()[0]
                    prodInt = spec.get_peaks()[1]
                    # Create a ms2Spectrum object
                    temp = ms2Spectrum(precsMz=precsMz, rt=rt, prodMz=prodMz, prodInt=prodInt)
                    self.ms2Data.append(temp)


    def findTopSignals(self, parameters):
        """
        Function to find the top signals in the raw LCMS data.

        Parameters
        ----------------------------------------------------------
        parameters: dict
            A dictionary containing the parameters for finding the top
            signals.
        """

        self.preFindTopSignals(parameters)
        self.purifyTopSignals(precsMzTol= parameters["precsMzTol"])
        self.deisotope()


    def preFindTopSignals(self, parameters):
        """
        Function to find metabolic features with the highest MS intensities.
        This function defines self.topSignals to store the selected.
        
        Parameters
        ----------------------------------------------------------
        parameters: dict
            A dictionary containing the parameters for finding the top
            signals.
        """

        # Set parameters
        sNum = parameters["sNum"]
        precsMzTol = parameters["precsMzTol"]
        prodMzTol = parameters["prodMzTol"]
        rtTol = parameters["rtTol"]
        intTol = parameters["intTol"]

        # Store the number of top signals
        self.sNum = sNum
        # Find top signals
        sNum = sNum * 2

        topMz = np.repeat(0.0, sNum)
        topInt = np.repeat(0.0, sNum)
        topRT = np.repeat(0.0, sNum)
        topScanNumber = np.repeat(0, sNum)
        for sidx, s in enumerate(self.ms1Data):
            mz = s.mz
            intensity = s.intensity
            rt = s.rt

            mz = mz[intensity > intTol]
            intensity = intensity[intensity > intTol]

            for idx, i in enumerate(intensity):
                m = mz[idx]
                isMatched, matchIdx = isFeatFound(topMz, topRT, m, rt, precsMzTol, rtTol*60)
                if isMatched:
                    if i > topInt[matchIdx]:
                        topMz[matchIdx] = m
                        topInt[matchIdx] = i
                        topRT[matchIdx] = rt
                        topScanNumber[matchIdx] = sidx
                else:
                    if i > np.amin(topInt):
                        newIdx = np.argmin(topInt)
                        topMz[newIdx] = m
                        topInt[newIdx] = i
                        topRT[newIdx] = rt
                        topScanNumber[newIdx] = sidx

        # Store the selected top signals as the MetaFeature objects
        self.topSignals = []
        for i in range(len(topMz)):
            f = MetaFeature(mz=topMz[i], rt=topRT[i],
                            intensity=topInt[i], scanIdx=topScanNumber[i])
            # Find the MS2 spectrum for top signals
            f.findMS2(ms2Data=self.ms2Data, mzTol=prodMzTol)
            self.topSignals.append(f)

        # Sort the top signals by mz
        self.topSignals = sorted(self.topSignals, key=lambda x: x.mz)


    def purifyTopSignals(self, precsMzTol=0.01, rtTol=0.5, fd=0.5):
        """
        A simple function to remove the top signals with poor peak shapes.

        Parameters
        ----------------------------------------------------------
        precsMzTol: float
            Tolerance of the precursor m/z difference.
        rtTol: float
            Retention time window (+- this value) for average calculation.
        fd: float
            The average intensity in retention time window needs to be lower 
            than fd * feature intensity.
        """

        # Find the top signals with good peak shapes
        purified = []
        p = len(self.topSignals)

        # Loop through all top signals
        for f in self.topSignals:
            # Find the MS1 scans within the retention time window
            ms1Range = [s for s in self.ms1Data if f.rt - rtTol < s.rt < f.rt + rtTol]
            intSeq = np.array([])
            mz = f.mz

            # Loop through all MS1 scans
            for s in ms1Range:
                mzScan = s.mz
                intensityScan = s.intensity
                mzdiff = abs(mzScan - mz)

                if any(mzdiff < precsMzTol):
                    intSeq = np.append(intSeq, intensityScan[np.argmin(mzdiff)])
                else:
                    intSeq = np.append(intSeq, 0)
            
            # Check the average intensity in retention time window
            if np.mean(intSeq) < fd * f.intensity:
                purified.append(f)
        
        self.topSignals = purified

        print("Top signal number reduced from {} to {} after purification.".format(p, len(self.topSignals)))


    def deisotope(self, precsMzTol=0.005, rtTol=0.001):
        """
        Function to remove the isotope from top signals.

        Parameters
        ----------------------------------------------------------
        precsMzTol: float
            Tolerance of the precursor m/z difference.
        rtTol: float
            Retention time window (+- this value) for average calculation.
        """

        # Store the number of top signals
        p = len(self.topSignals)

        # Check the number of isotopic peaks in top signals
        mzSeq = np.array([f.mz for f in self.topSignals])
        rtSeq = np.array([f.rt for f in self.topSignals])
        od = np.argsort(mzSeq)
        mzSeq = mzSeq[od]
        rtSeq = rtSeq[od]

        deIso = np.ones(len(self.topSignals), dtype=bool)
        for idx, _ in enumerate(mzSeq[:-1]):
            mzDiffBool = abs(mzSeq[idx + 1:] - mzSeq[idx] - 1.003355) < precsMzTol
            rtDiffBool = abs(rtSeq[idx + 1:] - rtSeq[idx]) < rtTol
            if any(mzDiffBool & rtDiffBool):
                t = idx + np.argmax(mzDiffBool & rtDiffBool) + 1
                deIso[t] = False

        self.topSignals = [i for idx, i in enumerate(self.topSignals) if deIso[idx] == True]
        print("Top signal number reduced from {} to {} after deisotoping.".format(p, len(self.topSignals)))

        # Sort the top signals by intensity from high to low
        self.topSignals = sorted(self.topSignals, key=lambda x: x.intensity, reverse=True)
        # Only keep a certain number of the highest intensity features
        self.topSignals = self.topSignals[:self.sNum]


    def computeSepEff(self, rtRange):
        """
        Function to compute the separation efficiency using top signals.

        Parameters
        ----------------------------------------------------------
        rtRange: list
            Retention time range for the separation efficiency calculation.
        """

        rtSeq = np.array([f.rt for f in self.topSignals])
        self.sepEff = sepEfficiency(rtSeq, rtRange)


    def getBPCData(self):
        """
        Function to prepare data for BPC.
        """
        self.rtBPC = []
        self.intensityBPC = []
        for m in self.ms1Data:
            self.rtBPC.append(m.rt)
            self.intensityBPC.append(np.max(m.intensity))
        self.rtBPC = np.array(self.rtBPC)
        self.intensityBPC = np.array(self.intensityBPC)
        self.intensityBPC = self.intensityBPC / np.max(self.intensityBPC) * 100


    def plotBPC(self, pltName):
        """
        Function to plot BPC.

        Parameters
        ----------------------------------------------------------
        pltName: str
            Name of the plot.
        """

        plt.plot(self.rtBPC, self.intensityBPC, color="black", linewidth=1)
        plt.xlabel("Retention time, min", fontname='Arial', fontsize=15, labelpad=5)
        plt.ylabel("Intensity, %", fontname='Arial', fontsize=15, labelpad=5)
        plt.xticks(fontname='Arial', fontsize=10)
        plt.yticks(fontname='Arial', fontsize=10)
        plt.savefig(fname=pltName, bbox_inches='tight')
        plt.close()


class MetaFeature:
    """
    Class represents the metabolomics feature in a single file.
    A feature has properties including:
        m/z, retention time, intensity, MS2.
    """

    def __init__(self, mz, rt, intensity, scanIdx):
        """
        Function to initiate a MetaFeature object.

        Parameters
        ----------------------------------------------------------
        mz: float
            m/z value.
        rt: float
            Retention time.
        intensity: float
            Peak intensity.
        scanIdx: int
            MS1 scan index.
        """

        self.mz = mz
        self.rt = rt
        self.intensity = intensity
        self.scanIdx = scanIdx
        self.ms2 = ms2Spectrum()

    def findMS2(self, ms2Data, rtTol=1.0, mzTol=0.02):
        """
        Function to find the MS2 spectrum and store in self.ms2 (ms2Spectrum object).

        Parameters
        ----------------------------------------------------------
        ms2Data: list
            A list of ms2Spectrum objects representing MS2 spectra.
        rtTol: float
            Retention window to find MS2 spectrum (+/- rtTol), in minute.
        mzTol: float
            m/z tolerance.
        """

        # Find the MS2 spectrum using retention time
        ms2Range = [s for s in ms2Data if (self.rt - rtTol) < s.rt < (self.rt + rtTol)]
        currentInt = 0
        targetIdx = -1
        for idx, s in enumerate(ms2Range):
            if abs(s.precsMz - self.mz) < mzTol:
                if np.sum(s.prodInt) > currentInt:
                    currentInt = np.sum(s.prodInt)
                    targetIdx = idx
        if targetIdx >= 0:
            self.ms2 = ms2Range[targetIdx]

    def addMS2(self, ms2):
        """
        Function to set the MS2 spectrum for a feature.

        Parameters
        ----------------------------------------------------------
        ms2: ms2Spectrum object
            The MS2 spectrum being added.
        """

        self.ms2 = ms2


    def showFeatInfo(self):
        """
        Function to print the feature information.
        """

        print("m/z: " + str(np.round(self.mz, decimals=4)))
        print("rt: " + str(np.round(self.rt, decimals=2)) + " min.")
        print("Intensity: " + str(np.round(self.intensity, decimals=0)))


class ms1Spectrum:
    """
    Class represents the MS1 spectrum.
    A MS1 spectrum has properties including:
        source file, retention time, 
        precursor m/z and intensities.
    """

    def __init__(self, mz, rt, intensity):
        """
        Function to initiate the msmsSpectrum by precursor mz,
        retention time, and source file.

        Parameters
        ----------------------------------------------------------
        mz: numpy array
            m/z values.
        rt: float
            Retention time.
        intensity: numpy array
            Intensities of ions.
        """

        self.mz = mz
        self.rt = rt
        self.intensity = intensity


class ms2Spectrum:
    """
    Class represents the MS2 spectrum.
    A MS2 spectrum has properties including:
        precursor m/z, source file, retention time,
        product m/z values, product intensities.
    """

    def __init__(self, precsMz=np.nan, rt=np.nan, prodMz=np.array([]), prodInt=np.array([])):
        """
        Function to initiate the msmsSpectrum by precursor mz,
        retention time, and source file.

        Parameters
        ----------------------------------------------------------
        precsMz: float
            Precursor m/z value.
        rt: float
            Retention time in minute.
        """

        self.precsMz = precsMz
        self.rt = rt
        self.prodMz = prodMz
        self.prodInt = prodInt


    def showSpecInfo(self, showMS2=False):
        """
        Function to print the MS2 spectrum information.

        Parameters
        ----------------------------------------------------------
        showMS2: bool
            Whether to show MS2 information.
        """

        print("Precursor m/z: " + str(np.round(self.precsMz, decimals=4)))
        print("rt: " + str(np.round(self.rt, decimals=2)) + " min.")
        if showMS2:
            print('prodMz: ' + str(np.round(self.prodMz, decimals=4)))
            print("Intensity: " + str(np.round(self.prodInt, decimals=0)))


def isFeatFound(mzSeq, rtSeq, mz, rt, mzTol=0.01, rtTol=1.0):
    """
    Function to determine if a new mz value has been found in 
    the selected features
    Return: [bool, int]
        bool: Ture for matched.
        int: Matching index. -1 for not matched.

    Parameters
    ----------------------------------------------------------
    mzSeq: numpy.array
        m/z values of the selected features.
    rtSeq: numpy.array
        Retention times of the selected features.
    mz: float
        m/z value of the feature to be matched.
    rt: float
        Retention time of the feature to be matched. In minute.
    mzTol: float
        Tolerance of the m/z difference.
    rtTol: float
        Tolerance of the retention time difference. In minute.
    """

    mzdiff = abs(mzSeq - mz) < mzTol  # Bool array
    rtdiff = abs(rtSeq - rt) < rtTol  # Bool array

    for i in range(len(mzdiff)):
        if mzdiff[i] and rtdiff[i]:
            return [True, i]

    return [False, -1]


def drawEICAround(ms1Data, mz, rt, outDir, pltName="TBD.png",
                  title="TBD", output=True, mzTol=0.01, rtTol=1):
    """
    Function to draw the EIC for given mz value and retention time.

    Parameters
    ----------------------------------------------------------
    ms1Data: list
        A list of ms1Spectrum object representing ms1 spectra.
    mz: float
        m/z value of the feature to be matched.
    rt: float
        Retention time of the feature to be matched.
    outDir: str
        Directory to output the EIC.
    pltName: str
        File name of the output figure.
    title: str
        Title of this figure.
    output: boolean
        Whether to output the EIC.
    mzTol: float
        Tolerance of the m/z difference.
    rtTol: float
        Tolerance of the retention time difference. In second.
    """

    ms1Range = [s for s in ms1Data if rt - rtTol < s.rt < rt + rtTol]
    intSeq = np.array([])
    rtSeq = np.array([])

    for s in ms1Range:
        mzScan = s.mz
        intensityScan = s.intensity
        rtScan = s.rt

        mzdiff = abs(mzScan - mz)

        if any(mzdiff < mzTol):
            intSeq = np.append(intSeq, intensityScan[np.argmin(mzdiff)])
            rtSeq = np.append(rtSeq, rtScan)
        else:
            intSeq = np.append(intSeq, 0)
            rtSeq = np.append(rtSeq, rtScan)

    intSeq = np.insert(np.array([0, 0]), 1, intSeq)

    diff = rtSeq[1] - rtSeq[0]
    rtSeq = np.concatenate(([rtSeq[0] - diff], rtSeq, [rtSeq[-1] + diff]))

    if output:
        os.chdir(outDir)
        plt.plot(rtSeq, intSeq)
        plt.xlabel('retention time, min')
        plt.ylabel('intensity')
        plt.title(title)
        plt.savefig(pltName)
        plt.close()
    else:
        return [rtSeq, intSeq]


def dotProd(ms2A, ms2B, mzTol=0.02, rep='simple'):
    """
    Function to calculate the MS2 similarity by dot product.
    Return 2.0 if MS2 does not exist.

    Parameters
    ----------------------------------------------------------
    ms2A: ms2Spectrum object
        The first MS2.
    ms2B: ms2Spectrum object
        The second MS2.
    mzTol: float
        Tolerance of the m/z difference.
    rep: str
        Representation of the MS2 spectrum. 'simple' for only
        reporting the dot product. 'full' for reporting the
        dot product and number of matched peaks.
    """

    # Check if MS2 exist in both spectra
    if len(ms2A.prodMz) == 0 or len(ms2B.prodMz) == 0:
        return 0.0

    mz1 = np.copy(ms2A.prodMz)
    mz2 = np.copy(ms2B.prodMz)
    int1 = np.copy(ms2A.prodInt)
    int2 = np.copy(ms2B.prodInt)

    # Align two mass spectra
    diff = []

    for mz in mz2:
        diff.append(np.absolute(mz1 - mz))

    diff = np.matrix(diff)
    minDiff = diff.min()
    intTable = np.array([int1, np.zeros(len(int1))])
    while minDiff < mzTol:
        # find the index of that match
        row = len(mz1)
        mc = divmod(np.argmin(diff), row)
        intTable[1, mc[1]] = int2[mc[0]]
        diff[:, mc[1]] = 1
        diff[mc[0]] = 1
        mz2[mc[0]] = -1
        minDiff = diff.min()

    unmatched = int2[mz2 != -1]
    matchNum = len(int2[mz2 == -1])
    if len(unmatched) != 0:
        t = np.row_stack((np.zeros(len(unmatched)), unmatched))
        intTable = np.column_stack((intTable, t))

    dp = np.inner(intTable[0], intTable[1]) / (
            np.inner(intTable[0], intTable[0]) * np.inner(intTable[1], intTable[1])) ** 0.5

    if rep == 'simple':
        return dp
    elif rep == 'full':
        return {'dotprod': dp, 'matchNumber': matchNum}


def getUniqueMS2(d, rtTol=1.0, precsMzTol=0.01, dpTol=0.95, returnNum=False):
    """
    Function to calculate the number of unique MS2.

    Parameters
    ----------------------------------------------------------
    d: MSData object
        MS data
    rtTol: float
        Tolerance of retention time, in min.
    precsMzTol: float
        Tolerance of precursor m/z, in Da.
    dpTol: float
        Tolerance of dot product.
    returnNum: boolean
        True to return the number of unique MS2;
        False to return the list of unique MS2.

    Returns
    ----------------------------------------------------------
    List:
        Unique MS2 spectra (when returnNum is False).
    """

    uniqueMS2 = []
    temp = [{
        "ms2": [d.ms2Data[0]],
        "precsMz": d.ms2Data[0].precsMz,
        "bestms2": d.ms2Data[0]
    }]

    for ms2 in d.ms2Data[1:]:
        mzSeq = np.array([t["precsMz"] for t in temp])
        matchBoo = abs(ms2.precsMz - mzSeq) < precsMzTol
        if np.any(matchBoo):
            matchIdxes = [i for i in range(len(temp)) if matchBoo[i]]
            toCompare = [temp[t]["bestms2"] for t in matchIdxes]
            dpSeq = [dotProd(ms2, m) for m in toCompare]
            if max(dpSeq) > dpTol:
                matchIdx = matchIdxes[np.argmax(dpSeq)]
                temp[matchIdx]["ms2"].append(ms2)
                temp[matchIdx]["precsMz"] = np.mean(np.array([t.precsMz for t in temp[matchIdx]["ms2"]]))
                sumInt = np.array([np.sum(t.prodInt) for t in temp[matchIdx]["ms2"]])
                temp[matchIdx]["bestms2"] = temp[matchIdx]["ms2"][np.argmax(sumInt)]
            else:
                temp.append({
                    "ms2": [ms2],
                    "precsMz": ms2.precsMz,
                    "bestms2": ms2
                })
        else:
            temp.append({
                "ms2": [ms2],
                "precsMz": ms2.precsMz,
                "bestms2": ms2
            })

        trans = []
        for i, m in enumerate(temp):
            if ms2.rt - m["ms2"][-1].rt > rtTol:
                trans.append(i)
        for i in list(reversed(trans)):
            uniqueMS2.append(temp.pop(i))

    uniqueMS2 += temp

    if returnNum:
        return len(uniqueMS2)
    else:
        return uniqueMS2


def getUniqueMz(d, precsMzTol=0.01):
    """
    Function to calculate the number of unique m/z
    ----------------------------------------------------------
    d: MSData object
        MS data
    precsMzTol: float
        Tolerance of precursor m/z, in Da.

    Returns
    ----------------------------------------------------------
    List:
        Unique MS2 spectra.
    """

    allMz = [ms1.mz for ms1 in d.ms1Data]
    # merge sublist in allMz to a list
    allMz = [item for sublist in allMz for item in sublist]
    allMz.sort()

    # get all values in allMz with difference between each two less than precsMzTol
    uniqueMz = [allMz[0]]
    for i in range(1, len(allMz)):
        if allMz[i] - uniqueMz[-1] > precsMzTol:
            uniqueMz.append(allMz[i])

    return uniqueMz


def numberChimericMS2(d, precsMzDiff=5.0, precsTol=0.01, intTol=2.0):
    """
    Function to calculate the number of chimeric MS2.

    Parameters
    ----------------------------------------------------------
    d: MSData object
        MS data
    precsMzDiff: float
        m/z window for precursor ion selection.
    precsTol: float
        m/z tolerance of precursor ion.
    intTol: float
        Intensity tolerance of ions with similar m/z, in %.

    Returns
    ----------------------------------------------------------
    int:
        number of chimeric MS2.
    """

    uniqueMS2 = getUniqueMz(d)
    chiMS2 = []
    for i in uniqueMS2:
        ms2 = i['bestms2']
        if len(ms2.prodMz) != 0:
            mz = ms2.prodMz[ms2.prodInt / np.max(ms2.prodInt) * 100 > intTol]
            temp1 = np.abs(ms2.precsMz - mz) < precsMzDiff
            temp2 = np.abs(ms2.precsMz - mz) > precsTol
            temp3 = np.abs(ms2.precsMz - mz + 1.003355) > precsTol
            temp4 = np.abs(ms2.precsMz - mz + 2 * 1.003355) > precsTol
            temp = np.logical_and(temp1, temp2)
            temp = np.logical_and(temp, temp3)
            temp = np.logical_and(temp, temp4)
            if np.sum(temp) != 0:
                chiMS2.append(ms2)

    return chiMS2


def getMobilePhasePct(gradient, timePoints):
    """
    Function to calculate the percentage of mobile phase A used 
    through the gradient, which can be used to estimate the
    speed of gradient.

    Parameters
    ----------------------------------------------------------
    gradient: numpy array
        Sequence of mobile phase percentages at all time points.
    timePoints:
        Sequence of time points.

    Returns
    ----------------------------------------------------------
    float:
        The percentage of mobile phase used.
    """

    return np.round(np.trapz(gradient, timePoints) / (100 * (timePoints[-1] - timePoints[0])) * 100, 2)


def outputConfig(name, timePoints, gradient):
    """
    Function to output the gradient to csv file.

    Parameters
    ----------------------------------------------------------
    name: str
        Name of this LC gradient.
    gradient:numpy array
        Sequence of mobile phase percentages at all time points.
    timePoints: numpy array
        Sequence of time points.
    """

    fname = name + '.csv'
    df = pd.DataFrame(data=[timePoints, gradient, 100 - gradient], index=['RT', 'MP', 'counterMP'])
    df = pd.DataFrame.transpose(df)
    df.to_csv(fname, index=False)


def outputGradientFig(name, gradient, timePoints):
    """
    Function to output the gradient figure in png.

    Parameters
    ----------------------------------------------------------
    name: str
        Name of this LC gradient.
    gradient:numpy array
        Sequence of mobile phase percentages at all time points.
    timePoints: numpy array
        Sequence of time points.
    """

    plt.plot(timePoints, gradient, color="black", linewidth=1)
    plt.xlabel("Retention time, min", fontname='Arial', fontsize=10, labelpad=2)
    plt.ylabel("Mobile phase, %", fontname='Arial', fontsize=10, labelpad=2)
    plt.xticks(fontname='Arial', fontsize=10)
    plt.yticks(np.arange(0, 120, 20), fontname='Arial', fontsize=10)
    plt.savefig(fname=name, bbox_inches='tight')
    plt.close()


def sepEfficiency(rtSeq, rtRange):
    """
    Calculate the separation efficiency using a series of retention times.

    Parameters
    ----------------------------------------------------------
    rtSeq: numpy array
        Retention times of the top features.
    rtRange: numpy array
        Times when gradient begins and ends. In minute.

    Returns
    ----------------------------------------------------------
    float:
        Separation efficiency.
    """

    rtRange = np.array(rtRange, dtype=float)
    temp = np.insert(rtRange, 1, np.sort(rtSeq))
    timeItvl = temp[1:] - temp[:-1]
    sqrti = np.sum(timeItvl * timeItvl)
    eff = ((1 / sqrti) * (rtRange[1] - rtRange[0]) ** 2 - 1) / len(rtSeq)
    return eff


def computeSecondGradient(parameters, model):
    """
    Calculate the second gradient to run.

    Parameters
    ----------------------------------------------------------
    parameters: dict
        Global parameters.
    model: gpModel class
        Gaussian process model.

    Returns
    ----------------------------------------------------------
    numpy array:
        The second gradient.
    """

    g1 = parameters['grads']['Init_1']
    pct1 = getMobilePhasePct(g1, parameters['timePoints'])
    temp = np.min(model.gradientPct)
    if pct1 > temp:
        temp = (pct1 - temp) * 0.5 + temp
        for idx, pct in enumerate(model.gradientPct):
            if pct - temp > 0:
                gradientReturn = np.copy(g1)
                gradientReturn[parameters['isChangable']] = model.gridX[idx]
                return gradientReturn
    else:
        print("Warning, the lower bound of mobile phase percentage range is higher than the initial gradient.")
        return np.NaN


def calGradPoints(parameters):
    """
    Calculate the gradient points.

    Parameters
    ----------------------------------------------------------
    parameters: dict
        Parameters.
    """

    r = parameters['gradRange']
    s = parameters['gradStep']
    gradPoints = np.arange(r[0], r[1], s)
    gradPoints = np.append(gradPoints, r[1])

    parameters['gradPoints'] = gradPoints


def calTotalMP(gradient, timePoints, flowRate):
    """
    Calculate the total mobile phase used.

    Parameters
    ----------------------------------------------------------
    grad: numpy array
        The gradient.
    timePoints: numpy array
        Sequence of time points.
    flowRate: float
        Flow rate of the LC system.
    
    Returns
    ----------------------------------------------------------
    list:
        Total mobile phase used [strong, weak].
    """

    pct = getMobilePhasePct(gradient, timePoints)
    totalMP = flowRate * (timePoints[-1] - timePoints[0])
    strongMP = pct * totalMP / 100
    weakMP = totalMP - strongMP
    return [strongMP, weakMP]