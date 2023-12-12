<img title="archs4py" alt="archs4py" src="https://user-images.githubusercontent.com/32603869/243101679-c5147d56-fce0-4498-9577-a300df7d6dce.png">

# archs4py - Official Python package to load and query ARCHS4 data

Official ARCHS4 companion package. This package is a wrapper for basic H5 commands performed on the ARCHS4 data files. Some of the data access is optimized for specific query strategies and should make this implementation faster than manually querying the data. The package supports automated file download, mutithreading, and some convenience functions such as data normalization.

ARCHS4py also supports the ARCHS4 alignment pipeline. When aligning FASTQ files using ARCHS4py gene and transcript counts will be compatible with the preprocessed ARCHS4 samples.

[Installation](#installation) | [Download H5 Files](#usage) | [List H5 Contents](#list-data-fields-in-h5) | [Extract Counts](#data-access) | [Extract Meta Data](#meta-data) | [Normalize Samples](#normalizing-data) | [FASTQ Alignment](#sequence-alignment) | [Versions](#list-versions)

## ARCHS4 data

ARCHS4 data is regularly updated to include publically available gene expression samples from RNA-seq. ARCHS4 processes the major platforms for human and mouse. As of 6/2023 ARCHS4 encompasses more than 1.5 million RNA-seq samples. All samples in ARCHS4 are homogeniously processed. ARCHS4 does currently not decern whether samples are bulk or single-cell and purely crawls GEO. Since samples are not always correctly annotated as single cell ARCHS4 uses a machine learning approach to predict single-cell samples and associated a singlecellprobability to each sample. Samples with a value larger than 0.5 can be removed from the queries if needed.

## Installation

The python package can be directly installed from this GitHub repository using the following command (pip or pip3 depending on system setup)

```
pip3 install archs4py
```

## Usage

### Download data file


The data is stored in large HDF5 files which first need to be downloaded. HDF5 stores matrix information in a compressed datastructure that allows efficient data access to slices of the data. There are separate files for `human` and `mouse` data. The supported files are `gene counts` and `transcript counts`. As of 6/2023 the files are larger than 30GB and depending on the network speed will take some time to download.

```python
import archs4py as a4

file_path = a4.download.counts("human", path="", version="latest")
```

## List data fields in H5

The H5 files contain data and meta data information. To list the contents of ARCHS4 H5 files use the built in `ls` function.

```python
import archs4py as a4

file = "human_gene_v2.2.h5"
a4.ls(file)
```

## Data access

archs4py supports several ways to load gene expression data. When querying ARCHS4 be aware that when loading too many samples the system might run out of memory. (e.g. the meta data search term is very broad). In most cases loading several thousand samples at the same time should be no problem. To find relevant samples there are 5 main functions in the `archs4py.data` module. A function to extract N random samples `archs4py.data.rand()`, a function to extract samples by index `archs4py.data.index()`, a function to extract samples based on meta data search `archs4py.data.meta()`, a function to extract samples based on a list of geo accessions `archs4py.data.samples()` and lastly a function to extract all samples belonging to a series `archs4.data.series()`.

<span id="#extract-counts"></span>

#### Extract a random set of samples

To extract a random gene expression matrix use the `archs4py.data.rand()` function. The function will return a pandas dataframe with samples as columns and genes as rows.

```python
import archs4py as a4

#path to file
file = "human_gene_v2.2.h5"

# extract 100 random samples and remove sinle cell data
rand_counts = a4.data.rand(file, 100, remove_sc=True)

```

#### Extract samples at specified index positions

Extract samples based on their index positions in the H5 file.

```python
import archs4py as a4

#path to file
file = "human_gene_v2.2.h5"

# get counts for samples at position [0,1,2,3,4]
pos_counts = a4.data.index(file, [0,1,2,3,4])

```

#### Extract samples matching search term in meta data

The ARCHS4 H5 file contains all meta data of samples. Using meta data search all matching samples can be extracted with the use of search terms. There is also a `archs4py.meta` module that will only return meta data. Meta data fields to be returned can be specified `meta_fields=["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"]`

```python
import archs4py as a4

#path to file
file = "human_gene_v2.2.h5"

# search and extract samples matching regex (ignores whitespaces)
meta_counts = a4.data.meta(file, "myoblast", remove_sc=True)

```

#### Extract samples in a list of GEO accession IDs

Samples can directly be downloaded by providing a list of GSM IDs. Samples not contained in ARCHS4 will be ignored.

```python
import archs4py as a4

#path to file
file = "human_gene_v2.2.h5"

#get sample counts
sample_counts = a4.data.samples(file, ["GSM1158284","GSM1482938","GSM1562817"])

```

#### Extract samples belonging to a GEO series

To download all samples of a GEO series for example `GSE64016` use the series function.

```python
import archs4py as a4

#path to file
file = "human_gene_v2.2.h5"

#get sample counts
series_counts = a4.data.series(file, "GSE64016")

```

## Meta data

<span id="#extract-meta"></span>

Additinally to the data module archs4py also supports the extraction of meta data. It supports similar endpoints to the `archs4.data` module. Meta data fields can be specified with: `meta_fields=["geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title"]`

```python
import archs4py as a4

#path to file
file = "human_gene_v2.2.h5"

# get sample meta data based on search term
meta_meta = a4.meta.meta(file, "myoblast", meta_fields=["characteristics_ch1", "source_name_ch1"])

# get sample meta data
sample_meta = a4.meta.samples(file, ["GSM1158284","GSM1482938","GSM1562817"])

# get series meta data
series_meta = a4.meta.series(file, "GSE64016")

# get all entries of a meta data field for all samples. In this example get all sample ids and gene symbols in H5 file
all_samples = a4.meta.field(file, "geo_accession")
all_symbols = a4.meta.field(file, "symbol")
```

## Normalizing data
<span id="#normalize"></span>
The package also supports simple normalization. Currently supported are quantile normalization, log2 + quantile normalization, and cpm. In the example below we load 100 random samples and apply log quantile.

```python
import archs4py as a4

file = "human_gene_v2.2.h5"
rand_counts = a4.data.rand(file, 100)

#normalize using log quantile (method options for now = ["log_quantile", "quantile", "cpm", "tmm"])
norm_exp = a4.normalize(rand_counts, method="log_quantile")

```

## Sequence alignment
<span id="#align"></span>
The `align` module contains a replication of the ARCHS4 alignment pipeline. When used on FASTQ files the resulting gene or transcript counts are compatible with the previously aligned samples in ARCHS4. The package is highly automated and only required a path to a FASTQ file or a folder containing multiple FASTQ files. All file dependencies will downloaded automatically and index will be built when needed.

### Align FASTQ file

Pass either a single or paired FASTQ file. This function can return transcript count, gene counts, or transcript level TPM data.

```python

import archs4py as a4

a4.align.load(["SRR14457464"], "data/example_1")

result = a4.align.fastq("human", "data/example_1/SRR14457464.fastq", return_type="gene", identifier="symbol")

```

The next example is a SRR file that extracts into a pair of paired end FASTQ files. They can be passed to ARCHS4py like this:

```python
import archs4py as a4

# the sample is paired-end and will result in two files (SRR15972519_1.fastq, SRR15972519_2.fastq)
a4.align.load(["SRR15972519"], "data/example_2")

result = a4.align.fastq("mouse", ["data/example_2/SRR15972519_1.fastq", "data/example_2/SRR15972519_2.fastq"], return_type="transcript")

```

### Align FASTQ files from folder

Align all FASTQ files in folder using the function `a4.align.folder()`. ARCHS4py will automatically matching samples if data is paired end.

```python

import archs4py as a4

a4.align.load(["SRR15972519", "SRR15972520", "SRR15972521"], "data/example_3")

result = a4.align.folder("mouse", "data/example_3", return_type="gene", identifier="symbol")

```

## List versions
<span id="#version"></span>
ARCHS4 has different versions to download from. Recommended is the default setting, which will download the latest data release.

```python
import archs4 as a4

print(a4.versions())

```

# Citation

When using ARCHS4 please cite the following reference:

Lachmann, Alexander, Denis Torre, Alexandra B. Keenan, Kathleen M. Jagodnik, Hoyjin J. Lee, Lily Wang, Moshe C. Silverstein, and Avi Ma’ayan. "Massive mining of publicly available RNA-seq data from human and mouse." Nature communications 9, no. 1 (2018): 1366.
https://www.nature.com/articles/s41467-018-03751-6




