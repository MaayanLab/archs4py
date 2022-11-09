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
import wget

import archs4py.utils

def gene_counts(species, path="", version="latest"):
    conf = archs4py.utils.get_config()
    try:
        fpath = wget.download(conf["GENE_COUNTS_"+species.upper()][version]["primary"], path)
        print("file downloaded to", fpath)
    except Exception:
        fpath = wget.download(conf["GENE_COUNTS_"+species.upper()][version]["fallback"], path)
        print("file downloaded to", fpath)
