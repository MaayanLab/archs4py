import numpy as np
import pandas as pd
from collections import Counter

import h5py as h5
import re
import numpy as np
import pandas as pd
import tqdm

def meta(file, search_term, meta_fields=["characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"], remove_sc=False, silent=False):
    """
    Search for samples in a file based on a search term in specified metadata fields.

    Args:
        file (str): The file path or object containing the data.
        search_term (str): The term to search for. Case-insensitive.
        meta_fields (list, optional): The list of metadata fields to search within.
            Defaults to ["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"].
        remove_sc (bool, optional): Whether to remove single-cell samples from the results.
            Defaults to False.
        silent (bool, optional): Print progress bar.

    Returns:
        pd.DataFrame: DataFrame containing the extracted metadata, with metadata fields as columns and samples as rows.
    """
    search_term = search_term.upper()
    with h5.File(file, "r") as f:
        meta = []
        idx = []
        mfields = []
        for field in tqdm.tqdm(meta_fields, disable=not silent):
            if field in f["meta"]["samples"].keys():
                try:
                    meta.append([x.decode("UTF-8").upper() for x in list(np.array(f["meta"]["samples"][field]))])
                    mfields.append(field)
                except Exception:
                    x=0
        meta = pd.DataFrame(meta, index=mfields ,columns=[x.decode("UTF-8").upper() for x in list(np.array(f["meta"]["samples"]["geo_accession"]))])
        for i in tqdm.tqdm(range(meta.shape[0]), disable=silent):
            idx.extend([i for i, item in enumerate(meta.iloc[i,:]) if re.search(search_term, item.upper())])
        if remove_sc:
            singleprob = np.where(np.array(f["meta/samples/singlecellprobability"]) < 0.5)[0]
            idx = sorted(list(set(idx).intersection(set(singleprob))))
        else:
            idx = sorted(list(set(idx)))
    return meta.iloc[:,idx].T

def field(file, field):
    gene_meta = []
    transcript_meta = []
    with h5.File(file, 'r') as f:
        sample_meta = list(f["meta/samples"])
    try:
        with h5.File(file, 'r') as f:
            gene_meta = list(f["meta/genes"])
    except Exception:
        x = 0
    try:
        with h5.File(file, 'r') as f:
            transcript_meta = list(f["meta/transcripts"])
    except Exception:
        x = 0
    with h5.File(file, 'r') as f:
        if field in sample_meta:
            return [x.decode("UTF-8").upper() for x in list(np.array(f["meta"]["samples"][field]))]
        elif field in gene_meta:
            return [x.decode("UTF-8").upper() for x in list(np.array(f["meta"]["genes"][field]))]
        elif field in transcript_meta:
            return [x.decode("UTF-8").upper() for x in list(np.array(f["meta"]["transcripts"][field]))]
        else:
            raise("specified field does not exist. Choose from supported sample meta fields or gene meta fields. List fields ysing archs4py.ls(filename) function")

def samples(file, samples, meta_fields=["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"], silent=False):
    """
    Extracts metadata for specified samples from an HDF5 file.

    Args:
        file (str): Path to the HDF5 file.
        samples (list): List of samples to extract metadata for.
        meta_fields (list, optional): List of metadata fields to extract. Defaults to ["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"].
        silent (bool, optional): If True, disables the progress bar. Defaults to False.

    Returns:
        pandas.DataFrame: DataFrame containing the extracted metadata, with metadata fields as columns and samples as rows.
    """
    samples = set(samples)
    with h5.File(file, "r") as f:
        meta = []
        mfields = []
        meta_samples = np.array([x.decode("UTF-8").upper() for x in list(np.array(f["meta"]["samples"]["geo_accession"]))])
        idx = [i for i,x in enumerate(meta_samples) if x in samples]
        for field in tqdm.tqdm(meta_fields, disable=not silent):
            if field in f["meta"]["samples"].keys():
                try:
                    meta.append([x.decode("UTF-8").upper() for x in list(np.array(f["meta"]["samples"][field][idx]))])
                    mfields.append(field)
                except Exception:
                    x=0
        meta = pd.DataFrame(meta, index=mfields ,columns=[x.decode("UTF-8").upper() for x in list(np.array(f["meta"]["samples"]["geo_accession"][idx]))])
        inter = meta.columns.intersection(set(samples))
    return meta.loc[:,inter].T

def series(file, series, meta_fields=["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"], silent=False):
    """
    Extracts metadata for specified series from an HDF5 file.

    Args:
        file (str): Path to the HDF5 file.
        series: Series to extract metadata for.
        meta_fields (list, optional): List of metadata fields to extract. Defaults to ["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"].
        silent (bool, optional): If True, disables the progress bar. Defaults to False.

    Returns:
        pandas.DataFrame: DataFrame containing the extracted metadata, with metadata fields as columns and samples as rows.
    """
    with h5.File(file, "r") as f:
        meta = []
        mfields = []
        meta_series = np.array([x.decode("UTF-8").upper() for x in list(np.array(f["meta"]["samples"]["series_id"]))])
        idx = [i for i,x in enumerate(meta_series) if x == series]
        for field in tqdm.tqdm(meta_fields, disable=not silent):
            if field in f["meta"]["samples"].keys():
                try:
                    meta.append([x.decode("UTF-8").upper() for x in list(np.array(f["meta"]["samples"][field][idx]))])
                    mfields.append(field)
                except Exception:
                    x=0
        meta = pd.DataFrame(meta, index=mfields ,columns=[x.decode("UTF-8").upper() for x in list(np.array(f["meta"]["samples"]["geo_accession"][idx]))])
    return meta.T

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
