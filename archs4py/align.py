import xalign
import archs4py
import biomart
import numpy as np
import pandas as pd

conf = archs4py.utils.get_config()
gene_mapping = {}
gene_mapping["homo_sapiens"] = None
gene_mapping["mus_musculus"] = None

def fastq(species, fastq, release="latest", t=8, overwrite=False, return_type="transcript", identifier="symbol"):
    if species == "mouse":
        species = "mus_musculus"
    elif species == "human":
        species = "homo_sapiens"
    result = xalign.align_fastq(species, fastq, release=conf["ALIGNMENT"][str(release)]["release"], t=t, noncoding=True, overwrite=overwrite)
    result.set_index("transcript", inplace=True)
    result.index = [x.split(".")[0] for x in result.index]
    if return_type == "gene":
        return aggregate(result.loc[:,"reads"], species, release, identifier)
    elif return_type == "tpm":
        return result.loc[:,"tpm"]
    else:
        return result.loc[:,"reads"]

def folder(species, folder, return_type="transcript", release="latest", overwrite=False, t=8, identifier="symbol"):
    if species == "mouse":
        species = "mus_musculus"
    elif species == "human":
        species = "homo_sapiens"
    
    gene_count, transcript_count = xalign.align_folder(species, folder, release=conf["ALIGNMENT"][str(release)]["release"], t=t, noncoding=True, overwrite=overwrite)
    del gene_count
    if return_type == "transcript":
        return transcript_count
    else:
        return aggregate(transcript_count, species, release, identifier)

def aggregate(transcript_count, species, release, identifier):
    if gene_mapping[species] is None:
        gene_mapping[species] = get_ensembl_mappings(species, release)
    trans = transcript_count.copy()
    trans.index = [x.split(".")[0] for x in transcript_count.index]
    trans.index = gene_mapping[species].loc[trans.index, "ensembl_gene"]
    trans = trans.groupby(trans.index).sum().astype(np.uint64)
    if identifier == "symbol":
        gm = gene_mapping[species].copy()
        gm.index = gm.loc[:, "ensembl_gene"]
        gm = gm[~gm.index.duplicated(keep='first')]
        trans.index = gm.loc[trans.index, "symbol"]
    return trans

def get_ensembl_mappings(species, release):  
    server = biomart.BiomartServer(conf["ALIGNMENT"][str(release)]["biomart"])
    if species == "mus_musculus":
        mart = server.datasets['mmusculus_gene_ensembl']
        attributes = ['ensembl_transcript_id', 'mgi_symbol', 'ensembl_gene_id', 'gene_biotype']
    else:
        mart = server.datasets['hsapiens_gene_ensembl']
        attributes = ['ensembl_transcript_id', 'hgnc_symbol', 'ensembl_gene_id', 'gene_biotype']                                                     
    # Get the mapping between the attributes                                    
    response = mart.search({'attributes': attributes})
    data = response.raw.data.decode('ascii')
    ensembl_ids = []
    # Store the data in a dict                                                  
    for line in data.splitlines():                                              
        line = line.split('\t')                                
        ensembl_ids.append(line)
    gene_map = pd.DataFrame(ensembl_ids)
    gene_map.index = gene_map.iloc[:,0]
    nn = np.where(gene_map.iloc[:,1] == "")[0]
    gene_map.iloc[nn, 1] = gene_map.iloc[nn, 2]
    gene_map.columns = ["ensembl_transcript", "symbol", "ensembl_gene", "biotype"]
    gene_map = gene_map[~gene_map.index.duplicated(keep='first')]
    return gene_map

def load(sras, outfolder):
    xalign.sra.load_sras(sras, outfolder)
    print("download complete")