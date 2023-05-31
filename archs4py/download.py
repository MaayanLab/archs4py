import wget
import sys
import requests
from tqdm import tqdm
import archs4py.utils

def bar_progress(current, total, width=80, update_interval=10):
    current_gb = current / (1024**3)  # Convert current bytes to GB
    total_gb = total / (1024**3)  # Convert total bytes to GB
    
    if current % (update_interval * 1024**2) == 0:  # Update progress every 10 MB
        progress_message = "Downloading: %d%% [%.2f GB / %.2f GB]" % (current / total * 100, current_gb, total_gb)
        sys.stdout.write("\r" + progress_message)
        sys.stdout.flush()

def gene_counts(species, path="", type="GENE_COUNTS", version="latest"):
    conf = archs4py.utils.get_config()

    try:
        fpath = wget.download(conf[type][species.upper()][version]["primary"], out=path, bar=bar_progress)
        print("file downloaded to", fpath)
    except Exception:
        fpath = wget.download(conf[type][species.upper()][version]["fallback"], out=path, bar=bar_progress)
        print("file downloaded to", fpath)

