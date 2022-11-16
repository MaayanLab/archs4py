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

def get_meta(file):
    f = h5.File(file, "r")
    meta = {}
    for field in list(f["meta"]["samples"].keys()):
        meta[field] = [x.decode("UTF-8") for x in list(np.array(f["meta"]["samples"][field]))]
    f.close()
    return meta

def get_meta_sample_field(file, field):
    f = h5.File(file, "r")
    meta_data = [x.decode("UTF-8") for x in list(np.array(f["meta"]["samples"][field]))]
    f.close()
    return meta_data

def get_meta_gene_field(file, field):
    f = h5.File(file, "r")
    meta_data = [x.decode("UTF-8") for x in list(np.array(f["meta"]["genes"][field]))]
    f.close()
    return meta_data