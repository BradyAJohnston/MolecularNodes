import gzip
import io
import os
from pathlib import Path
import requests

CACHE_DIR = Path(Path.home(), "MolecularNodesCache").expanduser()


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


class StructureDownloader:
    def __init__(self, cache: str | Path | None = CACHE_DIR):
        self.cache = str(Path(cache).absolute()) if cache else None
        if self.cache and not os.path.isdir(self.cache):
            os.makedirs(self.cache)

    def download(
        self,
        code: str,
        format: str = "cif",
        database: str = "rcsb",
    ) -> Path | io.BytesIO | io.StringIO:
        """Downloads a structure from the specified protein data bank in the given format.

        Parameters
        ----------
        code : str
            The code of the file to fetch. Supports both traditional 4-character codes
            and new format codes with 'pdb_' prefix (e.g., 'pdb_00009bdt').
        format : str
            The format of the file. Defaults to "cif".
            Must be one of ['cif', 'pdb', 'bcif'].
        database : str
            The database to fetch the file from.
            Defaults to 'rcsb'.

        Returns
        -------
        Union[Path, io.BytesIO, io.StringIO]
            If cache is enabled, returns the path to the cached file as a Path.
            If cache is disabled, returns either:
            - io.BytesIO for binary formats (bcif)
            - io.StringIO for text formats (cif, pdb)

        Raises
        ------
        ValueError
            If the specified format is not supported.
        FileDownloadPDBError
            If there is an error downloading the file from the database.
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
            file = os.path.join(self.cache, filename)
            if os.path.exists(file):
                return Path(file)
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
        "Get the URL for downloading the given file form a particular database."

        if database in ["rcsb", "pdb", "wwpdb"]:
            if format == "bcif":
                return f"https://models.rcsb.org/{code}.bcif"
            else:
                return f"https://files.rcsb.org/download/{code}.{format}"
        elif database == "alphafold":
            return self.get_alphafold_url(code, format)
        else:
            raise ValueError(f"Database {database} not currently supported.")

    def get_alphafold_url(self, code: str, format: str) -> str:
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
        """
        if format not in ["pdb", "cif", "bcif"]:
            raise ValueError(
                f"Format {format} not currently supported from AlphaFold databse."
            )

        url = f"https://alphafold.ebi.ac.uk/api/prediction/{code}"
        response = requests.get(url)
        data = response.json()[0]
        return data[f"{format}Url"]
