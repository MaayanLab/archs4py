# archs4py - Official Python package to load and query ARCHS4 data

Official ARCHS4 compagnion package. This package is a wrapper for basic H5 commands performed on the ARCHS4 data files. Some of the data access is optimized for specific query strategies and should make this implementation faster than manually querying the data. The package supports automated file download, mutithreading, and some convenience functions such as data normalization.

## Installation

The python package can be directly installed from this GitHub repository using the following command (pip or pip3 depending on system setup)

```
pip3 install git+https://github.com/MaayanLab/archs4py.git
```

## Usage

### Download data file

The data is stored in a large H5 file which first needs to be downloaded.

```python
import archs4 as a4

file_path = a4.download.gene_counts("human", path="", version="latest")
```

### Data access

archs4py supports several ways to load gene expression data. When querying ARCHS4 be aware that when loading too many samples the system might run out of memory. (e.g. the meta data search term is very broad). In most cases loading several thousand samples at the same time should be no problem.

```python
import archs4py as a4

#path to file
file = "human_gene_v2.2.h5"

# get counts for samples at position [0,1,2,3,4]
pos_counts = a4.data.index(file, [0,1,2,3,4])

# extract 100 random samples
# filterSingle=True will only retrieve bulk gene expression
rand_counts = a4.data.rand(file, 100, filterSingle=False)

# search and extract samples matching regex (ignores whitespaces)
# filterSingle=True will only retrieve bulk gene expression
meta_counts = a4.data.meta(file, "myoblast", filterSingle=False)

# get samples of specified series
series_counts = a4.data.series(file, "GSE150819")

#get sample counts
sample_counts = a4.data.samples(file, ["GSM1158284","GSM1482938","GSM1562817"])

```

### Direct access from S3

Gene expression can be loaded directly from S3 without downloading the complete file. This is very slow and will only work for few samples at a time. Instead of passing a file path pass the URL in S3.

```python
import archs4py as a4

#path to file
url = "https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.2.h5"

# extract 100 random samples
# filterSingle=True will only retrieve bulk gene expression
rand_counts = a4.data.rand(url, 100, filterSingle=False)

```

### Normalizing data

The package also supports simple normalization. Currently supported are quantile normalization, log2 + quantile normalization, and cpm

```python
import archs4py as a4

rand_counts = a4.data.rand(file, 100)

#normalize using log quantile (method options for now = ["log_quantile", "quantile", "cpm"])
norm_exp = a4.normalize(rand_counts, method="log_quantile")
```


### List versions

```python
import archs4 as a4

print(a4.versions())

```
