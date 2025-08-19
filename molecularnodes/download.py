import gzip
import io
from pathlib import Path
import requests

CACHE_DIR = Path(Path.home(), "MolecularNodesCache").expanduser()


class FileDownloadPDBError(Exception):
    """
    Exception raised for errors in the file download process.

    Attributes
    ----------
    message : str
        Explanation of the error
    """

    def __init__(
        self,
        message="There was an error downloading the file from the Protein Data Bank. PDB or format for PDB code may not be available.",
    ):
        self.message = message
        super().__init__(self.message)


class StructureDownloader:
    """
    A class for downloading molecular structure files from various databases.

    Parameters
    ----------
    cache : str, Path, or None, optional
        Directory path for caching downloaded files. If None, files are not cached.
        Defaults to CACHE_DIR.
    """

    def __init__(self, cache: str | Path | None = CACHE_DIR):
        if cache:
            cache_path = Path(cache).absolute()
            cache_path.mkdir(parents=True, exist_ok=True)
            self.cache = cache_path
        else:
            self.cache = None

    def download(
        self,
        code: str,
        format: str = "cif",
        database: str = "rcsb",
    ) -> Path | io.BytesIO | io.StringIO:
        """Download a structure from the specified protein data bank in the given format.

        Parameters
        ----------
        code : str
            The code of the file to fetch. Supports both traditional 4-character codes
            and new format codes with 'pdb_' prefix (e.g., 'pdb_00009bdt').
        format : str, optional
            The format of the file. Must be one of ['cif', 'pdb', 'bcif'].
            Defaults to "cif".
        database : str, optional
            The database to fetch the file from. Must be one of ['rcsb', 'pdb', 'wwpdb', 'alphafold'].
            Defaults to 'rcsb'.

        Returns
        -------
        Path or io.BytesIO or io.StringIO
            If cache is enabled, returns the path to the cached file as a Path.
            If cache is disabled, returns either:
            - io.BytesIO for binary formats (bcif)
            - io.StringIO for text formats (cif, pdb)

        Raises
        ------
        ValueError
            If the specified format is not supported, or if new format PDB codes
            are used with PDB format.
        FileDownloadPDBError
            If there is an error downloading the file from the database.

        Examples
        --------
        >>> downloader = StructureDownloader()
        >>> path = downloader.download("1abc", format="cif")
        >>> path = downloader.download("pdb_00009bdt", format="bcif")
        """
        code = code.strip()
        format = format.strip(".")
        supported_formats = ["cif", "pdb", "bcif"]
        if format not in supported_formats:
            raise ValueError(f"File format '{format}' not in: {supported_formats=}")

        # Check if the code has the new format prefix and is requesting PDB format
        if code.startswith("pdb_") and format == "pdb":
            raise ValueError(
                "New format PDB codes (starting with 'pdb_') are not compatible with .pdb format. Please use 'cif' or 'bcif' format instead."
            )

        _is_binary = format in ["bcif"]
        filename = f"{code}.{format}"

        if self.cache:
            file = self.cache / filename
            if file.exists():
                return file
        else:
            file = None

        try:
            r = requests.get(self._url(code, format, database))
            r.raise_for_status()
        except requests.HTTPError as e:
            raise FileDownloadPDBError(str(e))

        if _is_binary:
            content = r.content
            # Check if the content is gzipped
            if content[:2] == b"\x1f\x8b":  # gzip magic number
                content = gzip.decompress(content)
        else:
            content = r.text

        if file:
            mode = "wb+" if _is_binary else "w+"
            with open(file, mode) as f:
                f.write(content)
            return Path(file)
        else:
            if _is_binary:
                if not isinstance(content, bytes):
                    raise ValueError(
                        "Binary content is not bytes, please check your format."
                    )
                file = io.BytesIO(content)
            else:
                if not isinstance(content, str):
                    raise ValueError(
                        "Text content is not str, please check your format."
                    )
                file = io.StringIO(content)

        return file

    def _url(self, code: str, format: str, database: str = "rcsb") -> str:
        """Get the URL for downloading the given file from a particular database.

        Parameters
        ----------
        code : str
            The structure code to download.
        format : str
            The file format ('cif', 'pdb', or 'bcif').
        database : str, optional
            The database name. Defaults to 'rcsb'.

        Returns
        -------
        str
            The complete URL for downloading the file.

        Raises
        ------
        ValueError
            If the database is not currently supported.
        """

        if database in ["rcsb", "pdb", "wwpdb"]:
            if format == "bcif":
                return f"https://models.rcsb.org/{code}.bcif"
            else:
                return f"https://files.rcsb.org/download/{code}.{format}"
        elif database == "alphafold":
            return get_alphafold_url(code, format)
        else:
            raise ValueError(f"Database {database} not currently supported.")


def get_alphafold_url(code: str, format: str) -> str:
    """Get the URL for downloading a structure from AlphaFold database.

    Parameters
    ----------
    code : str
        The UniProt ID or AlphaFold DB identifier.
    format : str
        The file format to download ('pdb', 'cif', or 'bcif').

    Returns
    -------
    str
        The URL to download the structure file.

    Raises
    ------
    ValueError
        If the requested format is not supported.
    requests.RequestException
        If there is an error fetching data from the AlphaFold API.

    Examples
    --------
    >>> url = get_alphafold_url("P12345", "pdb")
    >>> print(url)
    https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.pdb
    """
    if format not in ["pdb", "cif", "bcif"]:
        raise ValueError(
            f"Format {format} not currently supported from AlphaFold database."
        )

    url = f"https://alphafold.ebi.ac.uk/api/prediction/{code}"
    response = requests.get(url)
    response.raise_for_status()  # This will raise an exception for HTTP errors
    data = response.json()[0]
    return data[f"{format}Url"]
