import os
import requests
import io

def fetch(code, format="cif", cache=None, database='rcsb'):
    """
    Downloads a structure from the specified protein data bank in the given format.

    Parameters
    ----------
    code : str
        The code of the file to fetch.
    format : str, optional
        The format of the file. Defaults to "cif". Possible values are ['cif', 'pdb', 
        'mmcif', 'pdbx', 'mmtf', 'bcif'].
    cache : str, optional
        The cache directory to store the fetched file. Defaults to None.
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
    supported_formats = ['cif', 'pdb', 'mmtf', 'bcif']
    if format not in supported_formats:
        raise ValueError(f"File format '{format}' not in: {supported_formats=}")

    _is_binary = (format in ['bcif', 'mmtf'])
    filename = f"{code}.{format}"
    # create the cache location
    if cache:
        if not os.path.isdir(cache):
            os.makedirs(cache)

        file = os.path.join(cache, filename)
    else:
        file = None

    # get the contents of the url
    r = requests.get(_url(code, format, database))
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
    
    if database == "rcsb":
        if format == "bcif":
            return f"https://models.rcsb.org/{code}.bcif"
        if format == "mmtf":
            return f"https://mmtf.rcsb.org/v1.0/full/{code}"
            
        else:
            return f"https://files.rcsb.org/download/{code}.{format}"
    # if database == "pdbe":
    #     return f"https://www.ebi.ac.uk/pdbe/entry-files/download/{filename}"
    else:
        ValueError(f"Database {database} not currently supported.")

