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

def gene_counts(species, version="latest"):
    conf = archs4py.utils.get_config()
    try:
        wget.download(conf["GENE_COUNTS_"+species.upper()][version]["primary"])
    except Exception:
        wget.download(conf["GENE_COUNTS_"+species.upper()][version]["fallback"])
