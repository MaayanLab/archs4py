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

def meta(file, search_term, meta_fields=["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"]):
    search_term = re.sub(r"_|-|'|/| |\.", "", search_term.upper())
    print("Searches for any occurrence of", search_term, "as regular expression")
    f = h5.File(file, "r")
    idx = []
    for field in meta_fields:
        if field in f["meta"]["samples"].keys():
            meta = [x.decode("UTF-8") for x in list(np.array(f["meta"]["samples"][field]))]
            idx.extend([i for i, item in enumerate(meta) if re.search(search_term, re.sub(r"_|-|'|/| |\.", "", item.upper()))])
    f.close()
    idx = sorted(list(set(idx)))
    counts = counts(file, idx)
    return counts

def rand(file, number, seed=1):
    random.seed(seed)
    f = h5.File(file, "r")
    gsm_ids = [x.decode("UTF-8") for x in np.array(f["meta/samples/geo_accession"])]
    f.close()
    idx = sorted(random.sample(range(len(gsm_ids)), number))
    return index(file, idx)

def series(file, series_id):
    f = h5.File(file, "r")
    series = [x.decode("UTF-8") for x in np.array(f["meta/samples/series_id"])]
    f.close()
    idx = [i for i,x in enumerate(series) if x == series_id]
    if len(idx) > 0:
        return index(file, idx)

def samples(file, sample_ids):
    sample_ids = set(sample_ids)
    f = h5.File(file, "r")
    samples = [x.decode("UTF-8") for x in np.array(f["meta/samples/geo_accession"])]
    f.close()
    idx = [i for i,x in enumerate(samples) if x in sample_ids]
    if len(idx) > 0:
        return index(file, idx)

def index(file, sample_idx, gene_idx = []):
    sample_idx = sorted(sample_idx)
    gene_idx = sorted(gene_idx)
    f = h5.File(file, "r")
    genes = np.array([x.decode("UTF-8") for x in np.array(f["meta/genes/gene_symbol"])])
    if len(sample_idx) == 0:
        return pd.DataFrame(index=genes[gene_idx])
    gsm_ids = np.array([x.decode("UTF-8") for x in np.array(f["meta/samples/geo_accession"])])[sample_idx]
    f.close()
    if len(gene_idx) == 0:
        gene_idx = list(range(len(genes)))
    exp = []
    PROCESSES = 16
    with multiprocessing.Pool(PROCESSES) as pool:
        results = [pool.apply_async(get_sample, (file, i, gene_idx)) for i in sample_idx]
        for r in tqdm.tqdm(results):
            res = r.get()
            exp.append(res)
    exp = np.array(exp).T
    exp = pd.DataFrame(exp, index=genes[gene_idx], columns=gsm_ids, dtype=np.uint32)
    return exp

def get_sample(file, i, gene_idx):
    try:
        f = h5.File(file, "r")
        temp = np.array(f["data/expression"][:,i], dtype=np.uint32)[gene_idx]
        f.close()
    except Exception:
        dd = np.array([0]*len(gene_idx))
        return dd
    return temp

