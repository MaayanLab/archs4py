import numpy as np
import pandas as pd
from collections import Counter

import h5py as h5
import tqdm
import re
import os
import json

import qnorm
import multiprocessing
import random

import archs4py.data
import archs4py.download
import archs4py.meta
import archs4py.utils

import importlib
importlib.reload(archs4py.data)
importlib.reload(archs4py.download)
importlib.reload(archs4py.meta)
importlib.reload(archs4py.utils)

def normalize(counts, method="log_quantiles"):
    norm_exp = 0
    if method == "quantiles":
        norm_exp = qnorm.quantile_normalize(np.array(counts))
    elif method == "log_quantiles":
        norm_exp = qnorm.quantile_normalize(np.log2(1+np.array(counts)))
    elif method == "cpm":
        g = np.array(counts)
        norm_exp = np.abs(g/g.sum(axis=0))*1_000_000
    norm_exp = pd.DataFrame(norm_exp, index=counts.index, columns=counts.columns, dtype=np.float32)
    return norm_exp

def get_meta(file):
    f = h5.File(file, "r")
    meta = {}
    for field in list(f["meta"]["samples"].keys()):
        meta[field] = [x.decode("UTF-8") for x in list(np.array(f["meta"]["samples"][field]))]
    f.close()
    return meta

def get_config():
    config_url = os.path.join(
        os.path.dirname(__file__),
        'data/config.json')
    with open(config_url) as json_file:
        data = json.load(json_file)
    return(data)

def versions():
    config = get_config()
    versions = config["GENE_COUNTS_HUMAN"].keys()
    return versions

