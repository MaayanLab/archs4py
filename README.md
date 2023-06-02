<img title="archs4py" alt="archs4py" src="https://github-production-user-asset-6210df.s3.amazonaws.com/32603869/242734021-a99ca725-6f10-4e01-85c1-3c1e1694dc68.png">

# archs4py - Official Python package to load and query ARCHS4 data

Official ARCHS4 compagnion package. This package is a wrapper for basic H5 commands performed on the ARCHS4 data files. Some of the data access is optimized for specific query strategies and should make this implementation faster than manually querying the data. The package supports automated file download, mutithreading, and some convenience functions such as data normalization.

## ARCHS4 data

ARCHS4 data is regularly updated to include publically available gene expression samples from RNA-seq. ARCHS4 processes the major platforms for human and mouse. As of 6/2023 ARCHS4 encompasses more than 1.5 million RNA-seq samples. All samples in ARCHS4 are homogeniously processed. ARCHS4 does currently not decern whether samples are bulk or single-cell and purely crawls GEO. Since samples are not always correctly annotated as single cell ARCHS4 uses a machine learning approach to predict single-cell samples and associated a singlecellprobability to each sample. Samples with a value larger than 0.5 can be removed from the queries if needed.

## Installation

The python package can be directly installed from this GitHub repository using the following command (pip or pip3 depending on system setup)

```
pip3 install git+https://github.com/MaayanLab/archs4py.git
```

## Usage

### Download data file

The data is stored in large HDF5 files which first need to be downloaded. HDF5 stores matrix information in a compressed datastructure that allows efficient data access to slices of the data. There are separate files for `human` and `mouse` data. The supported files are `gene counts` and `transcript counts`. As of 6/2023 the files are larger than 30GB and depending on the network speed will take some time to download.

```python
import archs4py as a4

file_path = a4.download.counts("human", path="", version="latest")
```

### Data access

archs4py supports several ways to load gene expression data. When querying ARCHS4 be aware that when loading too many samples the system might run out of memory. (e.g. the meta data search term is very broad). In most cases loading several thousand samples at the same time should be no problem. To find relevant samples there are 5 main functions in the `archs4py.data` module. A function to extract N random samples `archs4py.data.rand()`, a function to extract samples by index `archs4py.data.index()`, a function to extract samples based on meta data search `archs4py.data.meta()`, a function to extract samples based on a list of geo accessions `archs4py.data.samples()` and lastly a function to extract all samples belonging to a series `archs4.data.series()`.


#### Extract a random set of samples
```python
import archs4py as a4

#path to file
file = "human_gene_v2.2.h5"

# extract 100 random samples and remove sinle cell data
rand_counts = a4.data.rand(file, 100, remove_sc=True)
```

#### Extract samples at specified index positions
```python
import archs4py as a4

#path to file
file = "human_gene_v2.2.h5"

# get counts for samples at position [0,1,2,3,4]
pos_counts = a4.data.index(file, [0,1,2,3,4])

```

#### Extract samples matching search term in meta data
```python
import archs4py as a4

#path to file
file = "human_gene_v2.2.h5"

# search and extract samples matching regex (ignores whitespaces)
meta_counts = a4.data.meta(file, "myoblast", remove_sc=True)

```

#### Extract samples in a list of GEO accession ids
```python

#get sample counts
sample_counts = a4.data.samples(file, ["GSM1158284","GSM1482938","GSM1562817"])
```

#### Extract samples belonging to a GEO series
```python

#get sample counts
sample_counts = a4.data.series(file, ["GSM1158284","GSM1482938","GSM1562817"])
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
