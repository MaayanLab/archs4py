import wget
import sys
import requests
import archs4py.utils
import requests
import os

def bar_progress(current, total, width=80, update_interval=10):
    current_gb = current / (1024**3)  # Convert current bytes to GB
    total_gb = total / (1024**3)  # Convert total bytes to GB
    
    if current % (update_interval * 1024**2) == 0:  # Update progress every 10 MB
        progress_message = "Downloading: %d%% [%.2f GB / %.2f GB]" % (current / total * 100, current_gb, total_gb)
        sys.stdout.write("\r" + progress_message)
        sys.stdout.flush()

def counts(species, path="", type="GENE_COUNTS", version="latest"):
    """
    Download count files for a given species and count type.

    Args:
        species (str): The species for which count files are being downloaded. ["human", "mouse"]
        path (str, optional): The path where the downloaded file will be saved. Defaults to "".
        type (str, optional): The type of count file to be downloaded. Defaults to "GENE_COUNTS".
        version (str, optional): The version of the count file to be downloaded. Defaults to "latest". Versions can be listed with archs4py.versions()

    Returns:
        str: The path where the count file is downloaded.

    Raises:
        Exception: If an error occurs during the download process.

    Notes:
        The function first tries to download the count file using the primary URL specified in the configuration file.
        If the download fails, it falls back to the fallback URL specified in the configuration file.

        Supported count types:
        - GENE_COUNTS: Gene-level count files.
        - TRANSCRIPT_COUNTS: Transcript-level count files.
    """
    conf = archs4py.utils.get_config()

    try:
        file_name = os.path.basename(conf[type][species.upper()][version]["primary"])
        download_url = conf["DOWNLOAD_URL"]
        url = f"{download_url}?&file={file_name}&version=1337"
        response = requests.get(url)
    except:
        x = "just continue"
    
    try:
        fpath = wget.download(conf[type][species.upper()][version]["primary"], out=path, bar=bar_progress)
        print("file downloaded to", fpath)
    except Exception:
        fpath = wget.download(conf[type][species.upper()][version]["fallback"], out=path, bar=bar_progress)
        print("file downloaded to", fpath)

