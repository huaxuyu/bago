# imports
import numpy as np
import pandas as pd
import os

from .models import gpModel, compute_the_second_gradient
from .params import Parameters
from .raw_data_utils import evaluate_ms_data, output_gradient, read_method_from_csv


"""
project_path: str
    The path to the project directory. The directory should contain the following files:
    - data: folder containing the MS data
    - method: folder containing the LC gradient method files as .csv files

Example:

    my_project
    |- data
    |  |- 1.mzML
    |  |- 2.mzML
    |  |- ...
    |- method
    |  |- 1.csv
    |  |- 2.csv
    |  |- ...
"""

def initialize_model(p: Parameters):
    '''
    This function is used to initialize the model
    
    Parameters
    ----------

    Returns
    -------
    m: gpModel
        The model object
    p: Parameters
        The parameters object
    '''

    # STEP 1: file checking and parameter initialization
    # (1) check if the project path exists
    if not os.path.exists(p.project_dir):
        raise ValueError("The project path does not exist. Please check the path and try again.")
    # (2) check if the data folder exists
    if not os.path.exists(p.data_dir):
        raise ValueError("The data folder does not exist. Please create the 'data' folder for the raw MS data (mzML or mzXML).")
    # (3) check if the method folder exists
    if not os.path.exists(p.method_dir):
        raise ValueError("The method folder does not exist. Please create the 'method' folder for the LC gradient method files (.csv).")

    # STEP 2: read the initial raw MS data and run evaluation to start
    method_dic = _get_name_dic(p.method_dir)
    raw_names_dic = _get_name_dic(p.data_dir)
    if len(method_dic) == 0:
        raise ValueError("At least one method file is required to start the optimization.")
    elif len(raw_names_dic) == 0:
        raise ValueError("At least one raw MS data file is required to start the optimization.")
    elif len(method_dic) != len(raw_names_dic):
        raise ValueError("The number of method files and raw MS data files must be the same. Please check the files and try again.")
    else:
        for f in method_dic.keys():
            # (1) read the method file
            method_path = method_dic[f]
            p.time_arr, p.gradients[f] = read_method_from_csv(method_path)
            # (2) read the raw MS data file
            raw_data_path = raw_names_dic[f]
            exp = evaluate_ms_data(raw_data_path, p)
            exp.order = int(f)  # set the order of the experiment
            p.experiments.append(exp)

    # STEP 3: initialize the model
    m = gpModel()
    # Generate the search space
    print("Generating the search space...")
    m.generate_search_space(p)
    print("Search space contains {} different gradient programs.".format(len(m.gradient_arrs)))
    
    return m, p


def compute_next_gradient(m: gpModel, p: Parameters):
    '''
    Read new MS data and compute the next gradient.

    Parameters
    ----------
    m: gpModel
        The model object
    parameters: Parameters
        The parameters object
    '''

    # check if there is new MS data
    raw_names_dic = _get_name_dic(p.data_dir)
    
    exp_idx = [e.order for e in p.experiments]
    for k in raw_names_dic.keys():
        if int(k) not in exp_idx:
            raw_data_path = raw_names_dic[str(len(p.experiments) + 1)]
            exp = evaluate_ms_data(raw_data_path, p)
            p.experiments.append(exp)

    # if the second file is not provided, then just generate one using exploration
    if len(p.experiments) == 1:
        g = compute_the_second_gradient(parameters=p, model=m)
    else:
        # update the model and make predictions
        m.update_model(parameters=p)
        g = m.compute_next_gradient(acq=p.acq_func)
    
    name = str(len(p.experiments) + 1) + ".csv"
    path = os.path.join(p.method_dir, name)
    # output the gradient to a file
    output_gradient(path=path, time_arr=p.time_arr, pct_arr=g)
    p.gradients[str(len(p.experiments) + 1)] = g


# def runEvaluation(exp, parameters):
#     '''
#     This function is used to run the evaluation

#     Parameters
#     ----------------------------------------------------------
#     exp: dict
#         The experiment dictionary
#     parameters: dict
#         Global parameters.
#     '''
    
#     parameterInit(parameters)
#     readNewMSData(exp, parameters)

#     # Obtain the global separation index of each data
#     sepEffs = [exp[k].sepEff for k in exp.keys()]
#     parameters["sepEffs"] = {k: v for k, v in zip(exp.keys(), sepEffs)}

#     # Find the key of the data with the highest global separation index
#     maxSepEffKey = list(exp.keys())[np.argmax(sepEffs)]

#     # Find the number of unique MS/MS spectra
#     parameters['uniqueMS2'] = [raw_data_utils.getUniqueMS2(d=exp[k], returnNum=True) for k in exp.keys()]

#     if maxSepEffKey in parameters["grads"].keys():
#         # Find the gradient setting that gives the highest global separation index
#         maxSepEffGrad = parameters["grads"][maxSepEffKey]

#         # Print the gradient setting that gives the highest global separation index
#         print("The gradient setting that gives the highest global separation index is: {}.".format(maxSepEffGrad))

#     else:
#         print("The gradient setting that gives the highest global separation index is: {}.".format(maxSepEffKey))

#     # Create a data frame of two columns: gradient name, global separation index, and number of unique MS/MS spectra
#     df = pd.DataFrame({"Gradient": exp.keys(), "Global separation index": sepEffs, "Unique MS/MS": parameters['uniqueMS2']})

#     # Save the data frame to a csv file
#     df.to_csv("Evaluation.csv", index=False)


def _get_name_dic(path):
    raw_names = os.listdir(path)
    raw_names = [n for n in raw_names if not n.startswith(".")]
    raw_names_dic = {}
    for i in range(len(raw_names)):
        raw_names_dic[raw_names[i].split(".")[0]] = os.path.join(path, raw_names[i])
    return raw_names_dic