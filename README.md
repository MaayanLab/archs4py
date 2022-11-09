# archs4py - Python package to load and query ARCHS4 data

This package is a wrapper for basic H5 commands performed on the ARCHS4 data files. Some of the data access is optimized for specific query strategies and should make this implementation faster than manually querying the data. The package supports automated file download, mutithreading, and some convenience functions such as data normalization.

## Installation

The python package can be directly installed from this GitHub repository using the following command (pip or pip3 depending on system setup)

```
pip3 install git+https://github.com/MaayanLab/archs4py.git
```

## Usage

### Data access

archs4py supports several ways to load gene expression data. When querying ARCHS4 be aware that when loading too many samples the system might run out of memory. (e.g. the meta data search term is very broad). In most cases loading several thousand samples at the same time should be no problem.

```python
import archs4py as a4

#path to file
file = "archs4_gene_human_v2.1.2.h5"

# get counts for samples at position [0,1,2,3,4]
pos_counts = a4.get_counts(file, [0,1,2,3,4])

# extract 100 random samples
rand_counts = a4.get_random(file, 1000)

# search and extract samples matching regex (ignores whitespaces)
meta_counts = a4.meta_search(file, "myoblast")

# get samples of specified series
series_counts = a4.get_series(file, "GSE150819")

#normalize using log quantile (method options for now = ["log_quantile", "quantile", "cpm"])
norm_exp = a4.normalize(rand_counts, method="log_quantile")
```