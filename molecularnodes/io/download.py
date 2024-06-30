import io
import os
from pathlib import Path
import requests

CACHE_OLD = str(Path("~", ".MolecularNodes").expanduser())
CACHE_DIR = str(Path("~", "MolecularNodesCache").expanduser())

# rename old cache directories if users have them so we aren't leaving cached files in
# hidden folders on disk somewhere, I don't like the idea of silently renaming folders
# on a user's disk on load, so for now this will be disabled.
# TODO: make a decision on this (maybe a conformation popup on download)
# if os.path.exists(CACHE_OLD):
#     os.rename(CACHE_OLD, CACHE_DIR)


class FileDownloadPDBError(Exception):
    """
    Exception raised for errors in the file download process.

    Attributes:
        message -- explanation of the error
    """

    def __init__(
        self,
        message="There was an error downloading the file from the Protein Data Bank. PDB or format for PDB code may not be available.",
    ):
        self.message = message
        super().__init__(self.message)


def download(code, format="cif", cache=CACHE_DIR, database="rcsb"):
    """
    Downloads a structure from the specified protein data bank in the given format.

    Parameters
    ----------
    code : str
        The code of the file to fetch.
    format : str, optional
        The format of the file. Defaults to "cif". Possible values are ['cif', 'pdb',
        'mmcif', 'pdbx', 'bcif'].
    cache : str, optional
        The cache directory to store the fetched file. Defaults to `~/MolecularNodesCache`.
    database : str, optional
        The database to fetch the file from. Defaults to 'rcsb'.

    Returns
    -------
    file
        The fetched file as a file-like object.

    Raises
    ------
    ValueError
        If the specified format is not supported.
    """
    supported_formats = ["cif", "pdb", "bcif"]
    if format not in supported_formats:
        raise ValueError(f"File format '{format}' not in: {supported_formats=}")

    _is_binary = format in ["bcif"]
    filename = f"{code}.{format}"
    # create the cache location
    if cache:
        if not os.path.isdir(cache):
            os.makedirs(cache)

        file = os.path.join(cache, filename)
    else:
        file = None

    if file:
        if os.path.exists(file):
            return file

    # get the contents of the url
    try:
        r = requests.get(_url(code, format, database))
        r.raise_for_status()
    except requests.HTTPError:
        raise FileDownloadPDBError
    if _is_binary:
        content = r.content
    else:
        content = r.text

    if file:
        mode = "wb+" if _is_binary else "w+"
        with open(file, mode) as f:
            f.write(content)
    else:
        if _is_binary:
            file = io.BytesIO(content)
        else:
            file = io.StringIO(content)

    return file


def _url(code, format, database="rcsb"):
    "Get the URL for downloading the given file form a particular database."

    if database in ["rcsb", "pdb", "wwpdb"]:
        if format == "bcif":
            return f"https://models.rcsb.org/{code}.bcif"

        else:
            return f"https://files.rcsb.org/download/{code}.{format}"
    # if database == "pdbe":
    #     return f"https://www.ebi.ac.uk/pdbe/entry-files/download/{filename}"
    elif database == "alphafold":
        return get_alphafold_url(code, format)
    # if database == "pdbe":
    #     return f"https://www.ebi.ac.uk/pdbe/entry-files/download/{filename}"
    else:
        ValueError(f"Database {database} not currently supported.")


def get_alphafold_url(code, format):
    if format not in ["pdb", "cif", "bcif"]:
        ValueError(f"Format {format} not currently supported from AlphaFold databse.")

    # we have to first query the database, then they'll return some JSON with a list
    # of metadata, some items of which will be the URLs for the computed models
    # in different formats such as pdbUrl, cifUrl, bcifUrl
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{code}"
    response = requests.get(url)
    data = response.json()[0]
    return data[f"{format}Url"]
