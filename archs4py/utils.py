import numpy as np
import pandas as pd
from collections import Counter

import h5py as h5
import tqdm
import re
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

def normalize(counts, method="log_quantiles"):
    """
    Normalize the count matrix using a specified method.

    Args:
        counts (pd.DataFrame): A pandas DataFrame representing the count matrix.
        method (str, optional): The normalization method to be applied. Default is "log_quantiles".
            - "quantiles": Perform quantile normalization on the counts.
            - "log_quantiles": Perform quantile normalization on the log-transformed counts.
            - "cpm": Perform count per million (CPM) normalization.

    Returns:
        pd.DataFrame: A normalized count matrix as a pandas DataFrame with the same index and columns as the input.

    Raises:
        ValueError: If an unsupported normalization method is provided.
    """
    norm_exp = 0
    if method == "quantiles":
        norm_exp = qnorm.quantile_normalize(np.array(counts))
    elif method == "log_quantiles":
        norm_exp = qnorm.quantile_normalize(np.log2(1+np.array(counts)))
    elif method == "cpm":
        g = np.array(counts)
        norm_exp = np.abs(g/g.sum(axis=0))*1_000_000
    else:
        raise ValueError("Unsupported normalization method: " + method)
    norm_exp = pd.DataFrame(norm_exp, index=counts.index, columns=counts.columns, dtype=np.float32)
    return norm_exp