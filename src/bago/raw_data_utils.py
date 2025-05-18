# imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from masscube import process_single_file
from dataclasses import dataclass


@dataclass
class Experiment:
    """
    Class to store the experiment results.
    """

    order: int
    time_arr: np.ndarray
    pct_arr: np.ndarray
    total_mobile_phase_ratio: float
    total_feat_num: int
    good_feat_num: int
    is_good_feat: np.ndarray
    msms_num: int
    good_msms_num: int
    gaussian_sim_arr: np.ndarray
    noise_score_arr: np.ndarray
    pif_arr: np.ndarray
    bpc: np.ndarray


"""
Evaluation metrics:

1. Number of true peaks increased
2. Peak width
3. Peak quality by Gaussian similarity, noise score, and asymmetry factor
4. Number of unique MS/MS spectra on true peaks
5. MS/MS spectral quality by precursor ion fraction in MS1 scan
"""


def evaluate_ms_data(path, p):
    """
    Function to evaluate the MS data.

    Parameters
    ----------
    path: str
        Path to the MS data file.
    p: Parameters object
        Parameters.
    """

    d = process_single_file(path)

    order = len(p.experiments) + 1
    time_arr = p.time_arr
    pct_arr = p.gradients[str(order)]
    total_mobile_phase_ratio = cal_mobile_phase_ratio(time_arr, pct_arr)
    total_feat_num = len(d.features)
    is_good_feat = label_good_features(d, p.gaussian_sim_tol, p.noise_score_tol)
    good_feat_num = np.sum(is_good_feat)
    msms_num = len([f for f in d.features if f.ms2])
    good_msms_num = len([d.features[i] for i in range(len(d.features)) if d.features[i].ms2 and is_good_feat[i]])
    gaussian_sim_arr = np.array([f.gaussian_similarity for f in d.features])
    noise_score_arr = np.array([f.noise_score for f in d.features])
    pif_arr = np.array([cal_precursor_ion_fraction_using_MS1(d, f) for f in d.features])
    bpc = np.array([d.ms1_time_arr, d.base_peak_arr[:,1]])

    exp = Experiment(order=order,
                     time_arr=time_arr,
                     pct_arr=pct_arr,
                     total_mobile_phase_ratio=total_mobile_phase_ratio,
                     total_feat_num=total_feat_num,
                     good_feat_num=good_feat_num,
                     is_good_feat=is_good_feat,
                     msms_num=msms_num,
                     good_msms_num=good_msms_num,
                     gaussian_sim_arr=gaussian_sim_arr,
                     noise_score_arr=noise_score_arr,
                     pif_arr=pif_arr,
                     bpc=bpc)

    
    return exp


def label_good_features(d, gaussian_sim_tol=0.8, noise_score_tol=0.3):
    """
    Function to return a boolean array indicating whether the features are good or not.

    Parameters
    ----------
    d: dict
        Dictionary containing the MS data.
    gaussian_sim_tol: float
        Tolerance for Gaussian similarity.
    noise_score_tol: float
        Tolerance for noise score.

    Returns
    -------
    int:
        Number of good features.
    """

    a = np.zeros(len(d.features), dtype=bool)
    
    for i, feature in enumerate(d.features):
        if feature.gaussian_similarity > gaussian_sim_tol and feature.noise_score < noise_score_tol:
            a[i] = True
    
    return a


"""
Output and data visualization

=============================================

"""

def output_gradient(path, time_arr, pct_arr):
    """
    Function to output the gradient to csv file.

    Parameters
    ----------------------------------------------------------
    path: str
        Path to the csv file.
    time_arr: numpy array
        Time points where the mobile phase percentage is changed. Zero must be the first point. In minutes.
    pct_arr: numpy array
        Mobile phase percentage at each time point. In percent.
    """

    # c['Time', 'Target_mobile_phase', 'Complementary_mobile_phase']
    df = pd.DataFrame({
        'Time': time_arr,
        'Target_mobile_phase': pct_arr,
        'Complementary_mobile_phase': 100 - pct_arr
    })
    df.to_csv(path, index=False)


def outputGradientFig(name, time_arr, pct_arr, font_size=20, dpi=600, figsize=(6, 4)):
    """
    Function to output the gradient figure in png.

    Parameters
    ----------------------------------------------------------
    name: str
        Name of this LC gradient.
    time_arr: numpy array
        Time points where the mobile phase percentage is changed. Zero must be the first point. In minutes.
    pct_arr: numpy array
        Mobile phase percentage at each time point. In percent.
    font_size: int
        Font size of the figure.
    dpi: int
        Dots per inch of the figure.
    figsize: tuple
        Size of the figure.
    """

    if not name.endswith('.png'):
        name = name + '.png'

    plt.figure(figsize=figsize, dpi=dpi)
    plt.plot(time_arr, pct_arr, color="black", linewidth=1)
    plt.xlabel("Retention time, min", fontname='Arial', fontsize=font_size, labelpad=5)
    plt.ylabel("Mobile phase, %", fontname='Arial', fontsize=font_size, labelpad=5)
    plt.xticks(fontname='Arial', fontsize=font_size)
    plt.yticks(np.arange(0, 120, 20), fontname='Arial', fontsize=font_size)
    plt.savefig(fname=name, bbox_inches='tight')
    plt.close()


"""
Evaluation matrics for global compound separation
=============================================
"""


def gsi_from_rt_array(rt_array, rt_range):
    """
    Calculate the global separation index (GSI) from the retention time array.

    Parameters
    ----------
    rt_array: numpy array
        Retention time array of true peaks, in minutes.
    rt_range: list
        Retention time range as [start, end], in minutes.

    Returns
    -------
    float:
        The global separation index.
    """

    tmp = np.array([rt_range[0]] + list(rt_array) + [rt_range[1]])
    diff = np.diff(tmp)
    sqrti = np.sum(diff ** 2)   # sum of the squared retention time intervals
    t2 = (rt_range[1] - rt_range[0]) ** 2
    gsi = (t2/sqrti-1) / len(rt_array)
    
    return gsi


def overlap_index(rt_ranges):
    """
    Calculate the overlap index of retention time ranges, where 1 is no overlap
    and 0 is complete overlap. The overlap index belongs to the range (0, 1].
    
    Parameters
    ----------
    rt_ranges: list
        List of retention time ranges as [start, end], in minutes.

    Returns
    -------
    float:
        The overlap score.
    """

    rt_ranges = sorted(rt_ranges, key=lambda x: x[0])

    overlap = 0
    for i in range(len(rt_ranges)):
        for j in range(i + 1, len(rt_ranges)):
            if rt_ranges[i][1] > rt_ranges[j][0]:
                overlap += min(rt_ranges[i][1], rt_ranges[j][1]) - rt_ranges[j][0]
            else:
                break

    overlap = overlap / (rt_ranges[-1][1] - rt_ranges[0][0])

    # use exponential decay to penalize the overlap score
    overlap = np.exp(-overlap)
            
    return overlap


"""
Utility functions for computation
=============================================
"""


def cal_mobile_phase_ratio(time_arr, pct_arr):
    """
    Calculate the ratio of the mobile phase used during the run.

    Parameters
    ----------
    time_arr: numpy array
        Time points where the mobile phase percentage is changed. Zero must be the first point.
    pct_arr: numpy array
        Percentage of the mobile phase at each time point.

    Returns
    -------
    float:
        The ratio of mobile phase used.
    """

    return np.trapz(pct_arr, time_arr) / 100 / (time_arr[-1] - time_arr[0])


def cal_mobile_phase_volumes(time_arr, pct_arr, flow_rate):
    """
    Calculate the total mobile phase used.

    Parameters
    ----------------------------------------------------------
    time_arr: numpy array
        Time points where the mobile phase percentage is changed. Zero must be the first point. In minutes.
    pct_arr: numpy array
        Mobile phase percentage at each time point. 
    flow_rate: float
        Flow rate of the mobile phase, in mL/min.
    
    Returns
    ----------------------------------------------------------
    list:
        Total mobile phase used [strong, weak].
    """

    ratio = cal_mobile_phase_ratio(time_arr, pct_arr)
    total_volume = flow_rate * (time_arr[-1] - time_arr[0]) / 60

    return total_volume * ratio


def cal_precursor_ion_fraction_using_MS1(d, feature, isolation_window=1.0):
    """
    Calculate the precursor ion fraction (PIF) of a feature using MS1 data.

    Parameters
    ----------
    d: MSData object defined in masscube package
        MS data object containing the MS data.
    feature: Feature
        Feature object containing the MS data.
    isolation_window: float
        Isolation window for the precursor ion, in Da.

    Returns
    -------
    float:
        The precursor ion fraction.
    """

    if feature.ms2 is not None:
        mz = feature.ms2.precursor_mz
        rt = feature.ms2.time

        scan = d.find_ms1_scan_by_rt(rt)

        tmp = scan.signals
        tmp = tmp[np.abs(tmp[:, 0] - mz) < isolation_window, :]
        if len(tmp) > 0:
            k = np.argmin(np.abs(tmp[:, 0] - mz))
            return tmp[k, 1] / np.sum(tmp[:, 1])
        else:
            return 1.0
    else:
        return np.nan
    

def cal_precursor_ion_fraction_using_MS2(feature, isolation_window=1.0):
    """
    Calculate the precursor ion fraction (PIF) of a feature using the MS/MS spectrum.
    The MS/MS spectrum should not be cleaned (i.e. not removing ions > precursor ion m/z).

    Parameters
    ----------
    feature: Feature
        Feature object containing the MS data.
    isolation_window: float
        Isolation window for the precursor ion, in Da.

    Returns
    -------
    float:
        The precursor ion fraction.
    """

    if feature.ms2 is not None:
        tmp = feature.ms2.signals
        tmp = tmp[np.abs(tmp[:, 0] - feature.mz) < isolation_window, :]
        if len(tmp) > 0:
            k = np.argmin(np.abs(tmp[:, 0] - feature.mz))
            return tmp[k, 1] / np.sum(tmp[:, 1])
        else:
            return 1.0
    else:
        return np.nan


def read_method_from_csv(path):
    """
    Read the LC gradient method from a csv file.

    Parameters
    ----------
    path: str
        Path to the csv file.
    
    Returns
    -------
    time_arr: numpy array
        Time points where the mobile phase percentage is changed. Zero must be the first point. In minutes.
    pct_arr: numpy array
        Mobile phase percentage at each time point. In percent.
    """

    df = pd.read_csv(path)
    time_arr = df.iloc[:, 0].values
    pct_arr = df.iloc[:, 1].values

    return time_arr, pct_arr