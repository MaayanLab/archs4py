import numpy as np
import pandas as pd

import h5py as h5
import s3fs
import tqdm
import re

import multiprocessing
import random

def resolve_url(url):
    u1 = url.rsplit('/', 1)
    u2 = u1[0].rsplit('/', 1)
    file_name = u1[-1]
    bucket_name = u2[-1]
    endpoint = u2[0]
    S3_URL = "s3://"+bucket_name+"/"+file_name
    return(S3_URL, endpoint)

def fetch_meta_remote(field, s3_url, endpoint):
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': endpoint})
    with h5.File(s3.open(s3_url, 'rb'), 'r', lib_version='latest') as f:
        meta = [x.decode("UTF-8") for x in list(np.array(f[field]))]
    return np.array(meta)

def meta(file, search_term, meta_fields=["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"],  remove_sc=False, silent=False):
    """
    Search for samples in a file based on a search term in specified metadata fields.

    Args:
        file (str): The file path or object containing the data.
        search_term (str): The term to search for. The search is case-insensitive and supports regular expressions.
        meta_fields (list, optional): The list of metadata fields to search within.
            Defaults to ["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"].
        remove_sc (bool, optional): Whether to filter single-cell samples from the results.
            Defaults to False.
        silent (bool, optional): Print progress bar.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the gene expression data for the matching samples.
    """
    search_term = re.sub(r"_|-|'|/| |\.", "", search_term.upper())
    if not silent:
        print("Searches for any occurrence of", search_term, "as regular expression")
    if file.startswith("http"):
        return meta_remote(file, search_term, meta_fields, remove_sc, silent)
    else:
        return meta_local(file, search_term, meta_fields, remove_sc, silent)

def meta_local(file, search_term, meta_fields=["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"], remove_sc=False, silent=False):
    f = h5.File(file, "r")
    idx = []
    for field in meta_fields:
        if field in f["meta"]["samples"].keys():
            meta = [x.decode("UTF-8") for x in list(np.array(f["meta"]["samples"][field]))]
            idx.extend([i for i, item in enumerate(meta) if re.search(search_term, re.sub(r"_|-|'|/| |\.", "", item.upper()))])
    if remove_sc:
        singleprob = np.where(np.array(f["meta/samples/singlecellprobability"]) < 0.5)[0]
    f.close()
    if remove_sc:
        idx = sorted(list(set(idx).intersection(set(singleprob))))
    else:
        idx = sorted(list(set(idx)))
    counts = index(file, idx, silent=silent)
    return counts

def meta_remote(url, search_term, meta_fields=["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"], remove_sc=False, silent=False):
    s3_url, endpoint = resolve_url(url)
    idx = []
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': endpoint})
    with h5.File(s3.open(s3_url, 'rb'), 'r', lib_version='latest') as f:
        for field in meta_fields:
            if field in f["meta"]["samples"].keys():
                meta = [x.decode("UTF-8") for x in list(np.array(f["meta"]["samples"][field]))]
                idx.extend([i for i, item in enumerate(meta) if re.search(search_term, re.sub(r"_|-|'|/| |\.", "", item.upper()))])
        if remove_sc:
            singleprob = np.where(np.array(f["meta/samples/singlecellprobability"]) < 0.5)[0]
    if remove_sc:
        idx = sorted(list(set(idx).intersection(set(singleprob))))
    else:
        idx = sorted(list(set(idx)))
    counts = index_remote(url, idx, silent)
    return counts

def rand(file, number, seed=1, remove_sc=False, silent=False):
    """
    Randomly select a specified number of samples from a file.

    Args:
        file (str): The file path or object containing the data.
        number (int): The number of samples to select randomly.
        seed (int, optional): The seed value for the random number generator. Defaults to 1.
        remove_sc (bool, optional): Whether to remove single-cell samples from the selection. Defaults to False.
        silent (bool, optional): Print progress bar.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the randomly selected samples' gene expression data.
    """
    random.seed(seed)
    if file.startswith("http"):
        return rand_remote(file, number, remove_sc, silent)
    else:
        return rand_local(file, number, remove_sc, silent)

def rand_local(file, number, remove_sc, silent=False):
    f = h5.File(file, "r")
    gsm_ids = [x.decode("UTF-8") for x in np.array(f["meta/samples/geo_accession"])]
    if remove_sc:
        singleprob = np.array(f["meta/samples/singlecellprobability"])
    f.close()
    if remove_sc:
        idx = sorted(random.sample(list(np.where(singleprob < 0.5)[0]), number))
    else:
        idx = sorted(random.sample(range(len(gsm_ids)), number))
    return index(file, idx, silent=silent)

def rand_remote(url, number, remove_sc, silent=False):
    s3_url, endpoint = resolve_url(url)
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': endpoint})
    with h5.File(s3.open(s3_url, 'rb'), 'r', lib_version='latest') as f:
        number_samples = len(f["meta/samples/geo_accession"])
        if remove_sc:
            singleprob = np.array(f["meta/samples/singlecellprobability"])
    if remove_sc:
        idx = sorted(random.sample(list(np.where(singleprob < 0.5)[0]), number))
    else:
        idx = sorted(random.sample(range(number_samples), number))
    return index_remote(url, idx, silent=silent)

def series(file, series_id, silent=False):
    """
    Retrieve samples belonging to a specific series from a file.

    Args:
        file (str): The file path or object containing the data.
        series_id (str): The ID of the series to retrieve samples from.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the gene expression data for the samples belonging to the specified series.
    """
    if file.startswith("http"):
        return series_remote(file, series_id, silent=silent)
    else:
        return series_local(file, series_id, silent=silent)

def series_local(file, series_id, silent=False):
    f = h5.File(file, "r")
    series = [x.decode("UTF-8") for x in np.array(f["meta/samples/series_id"])]
    f.close()
    idx = [i for i,x in enumerate(series) if x == series_id]
    if len(idx) > 0:
        return index(file, idx, silent=silent)

def series_remote(url, series_id, silent=False):
    s3_url, endpoint = resolve_url(url)
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': endpoint})
    with h5.File(s3.open(s3_url, 'rb'), 'r', lib_version='latest') as f:
        series = [x.decode("UTF-8") for x in np.array(f["meta/samples/series_id"])]
    idx = [i for i,x in enumerate(series) if x == series_id]
    if len(idx) > 0:
        return index_remote(url, idx, silent)

def samples(file, sample_ids, silent=False):
    if file.startswith("http"):
        return samples_remote(file, sample_ids, silent=silent)
    else:
        return samples_local(file, sample_ids, silent=silent)

def samples_local(file, sample_ids, silent=False):
    sample_ids = set(sample_ids)
    f = h5.File(file, "r")
    samples = [x.decode("UTF-8") for x in np.array(f["meta/samples/geo_accession"])]
    f.close()
    idx = [i for i,x in enumerate(samples) if x in sample_ids]
    if len(idx) > 0:
        return index(file, idx, silent=silent)

def samples_remote(url, sample_ids, silent=False):
    sample_ids = set(sample_ids)
    s3_url, endpoint = resolve_url(url)
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': endpoint})
    with h5.File(s3.open(s3_url, 'rb'), 'r', lib_version='latest') as f:
        sample_ids = [x.decode("UTF-8") for x in np.array(f["meta/samples/geo_accession"])]
    idx = [i for i,x in enumerate(samples) if x in sample_ids]
    if len(idx) > 0:
        return index_remote(url, idx, silent=silent)

def index(file, sample_idx, gene_idx = [], silent=False):
    """
    Retrieve gene expression data from a specified file for the given sample and gene indices.

    Args:
        file (str): The file path or object containing the data.
        sample_idx (list): A list of sample indices to retrieve expression data for.
        gene_idx (list, optional): A list of gene indices to retrieve expression data for. Defaults to an empty list (return all).
        silent (bool, optional): Whether to disable progress bar. Defaults to False.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the gene expression data.
    """
    sample_idx = sorted(sample_idx)
    gene_idx = sorted(gene_idx)
    row_encoding = get_encoding(file)
    f = h5.File(file, "r")
    genes = np.array([x.decode("UTF-8") for x in np.array(f[row_encoding])])
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
        for r in tqdm.tqdm(results, disable=silent):
            res = r.get()
            exp.append(res)
    exp = np.array(exp).T
    exp = pd.DataFrame(exp, index=genes[gene_idx], columns=gsm_ids, dtype=np.uint32)
    return exp

def index_remote(url, sample_idx, gene_idx = [], silent=False):
    if len(sample_idx) == 0:
        return pd.DataFrame(index=genes[gene_idx])
    s3_url, endpoint = resolve_url(url)
    sample_idx = sorted(sample_idx)
    gene_idx = sorted(gene_idx)
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': endpoint})
    row_encoding = get_encoding_remote(s3, url)
    genes = fetch_meta_remote(row_encoding, s3_url, endpoint)
    if len(gene_idx) == 0:
        gene_idx = np.array(list(range(len(genes))))
    gsm_ids = fetch_meta_remote("meta/samples/geo_accession", s3_url, endpoint)[sample_idx]
    with h5.File(s3.open(s3_url, 'rb'), 'r', lib_version='latest') as f:
        exp = np.array(f["data/expression"][:,np.array(sample_idx)], dtype=np.uint32)[gene_idx]
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

def get_sample_remote(s3_url, endpoint, i, gene_idx):
    try:
        s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': endpoint})
        with h5.File(s3.open(s3_url, 'rb'), 'r', lib_version='latest') as f:
            temp = np.array(f["data/expression"][:,i], dtype=np.uint32)[gene_idx]
        return temp
    except Exception:
        dd = np.array([0]*len(gene_idx))
        return dd

def get_encoding(file):
    with h5.File(file) as f:
        if "genes" in list(f["meta"].keys()):
            if "gene_symbol" in list(f["meta/genes"].keys()):
                return "meta/genes/gene_symbol"
            elif "symbol" in list(f["meta/genes"].keys()):
                return "meta/genes/symbol"
        elif "transcripts" in list(f["meta"].keys()):
            if "ensembl_id" in list(f["meta/trancripts"].keys()):
                return "meta/trancripts/ensembl_id"
        else:
            raise Exception("error in gene/transcript meta data")

def get_encoding_remote(s3, s3_url):
    with h5.File(s3.open(s3_url, 'rb'), 'r', lib_version='latest') as f:
        if "genes" in list(f["meta"].keys()):
            if "gene_symbol" in list(f["meta/genes"].keys()):
                return "meta/genes/gene_symbol"
            elif "symbol" in list(f["meta/genes"].keys()):
                return "meta/gene_symbol"
        elif "transcripts" in list(f["meta"].keys()):
            if "ensembl_id" in list(f["meta/trancripts"].keys()):
                return "meta/trancripts/ensembl_id"
        else:
            raise Exception("error in gene/transcript meta data")