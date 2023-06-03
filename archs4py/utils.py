import numpy as np
import pandas as pd
from scipy import stats

import os
import json

import qnorm

def get_config():
    config_url = os.path.join(
        os.path.dirname(__file__),
        'data/config.json')
    with open(config_url) as json_file:
        data = json.load(json_file)
    return(data)


def versions():
    """
    Get the available versions of human gene counts.

    Returns:
        list: A list of available versions of human gene counts.
    """
    config = get_config()
    versions = config["GENE_COUNTS"]["HUMAN"].keys()
    return versions

def normalize(counts, method="log_quantile", tmm_outlier=0.05):
    """
    Normalize the count matrix using a specified method.

    Args:
        counts (pd.DataFrame): A pandas DataFrame representing the count matrix.
        method (str, optional): The normalization method to be applied. Default is "log_quantiles".
            - "quantile": Perform quantile normalization on the counts.
            - "log_quantile": Perform quantile normalization on the log-transformed counts.
            - "cpm": Perform count per million (CPM) normalization.
            - "tmm": Perform trimmed mean normalization

    Returns:
        pd.DataFrame: A normalized count matrix as a pandas DataFrame with the same index and columns as the input.

    Raises:
        ValueError: If an unsupported normalization method is provided.
    """
    norm_exp = 0
    if method == "quantile":
        norm_exp = qnorm.quantile_normalize(np.array(counts))
    elif method == "log_quantile":
        norm_exp = qnorm.quantile_normalize(np.log2(1+np.array(counts)))
    elif method == "cpm":
        g = np.array(counts)
        norm_exp = np.abs(g/g.sum(axis=0))*1_000_000
    elif method == "tmm":
        norm_exp = tmm_norm(counts, tmm_outlier)
    else:
        raise ValueError("Unsupported normalization method: " + method)
    norm_exp = pd.DataFrame(norm_exp, index=counts.index, columns=counts.columns, dtype=np.float32)
    return norm_exp

def tmm_norm(exp, percentage=0.05):
    lexp = np.log2(1+exp).astype(np.float32)
    tmm = trimmed_mean(lexp, percentage)
    nf = pd.DataFrame(np.tile(tmm, (exp.shape[0], 1)), index=lexp.index, columns=lexp.columns)
    temp = (lexp/nf)
    return temp

def trimmed_mean(matrix, percentage):
    matrix = np.array(matrix)
    trimmed_means = []
    for col in range(matrix.shape[1]):
        data = matrix[:, col].copy()
        data = data[data > 0]
        n_trim = int(len(data) * percentage)
        sorted_values = np.sort(data)
        trimmed_values = sorted_values[n_trim:-n_trim]
        trimmed_mean = np.mean(trimmed_values)
        trimmed_means.append(trimmed_mean)
    return trimmed_means