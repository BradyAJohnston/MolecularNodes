import gzip
import io
import os
from pathlib import Path
from typing import Literal, cast
from urllib.parse import quote
import requests
from biotite.database import afdb

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
            Defaults to "cif". Ignored for 'pubchem', which always downloads SDF.
        database : str, optional
            The database to fetch the file from. Must be one of ['rcsb', 'pdb', 'wwpdb', 'alphafold', 'pubchem'].
            Defaults to 'rcsb'. For 'pubchem' the code is a compound CID (e.g. '2244')
            or a compound name (e.g. 'aspirin'), and the 3D conformer is downloaded.

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
        if database == "pubchem":
            return self._download_pubchem(code)
        supported_formats = ["cif", "pdb", "bcif"]
        if format not in supported_formats:
            raise ValueError(f"File format '{format}' not in: {supported_formats=}")

        # Check if the code has the new format prefix and is requesting PDB format
        if code.startswith("pdb_") and format == "pdb":
            raise ValueError(
                "New format PDB codes (starting with 'pdb_') are not compatible with .pdb format. Please use 'cif' or 'bcif' format instead."
            )

        # Use biotite's afdb.fetch for AlphaFold database
        if database == "alphafold":
            # checked against `supported_formats` above
            afdb_format = cast(Literal["cif", "pdb", "bcif"], format)
            try:
                # biotite returns StringIO/BytesIO without a cache, or the file path
                # as a string with one
                if self.cache is None:
                    return afdb.fetch(code, format=afdb_format, overwrite=False)
                return Path(
                    afdb.fetch(
                        code,
                        format=afdb_format,
                        target_path=self.cache,
                        overwrite=False,
                    )
                )
            except Exception as e:
                raise FileDownloadPDBError(
                    f"Error downloading from AlphaFold database: {str(e)}"
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
            r = requests.get(self._url(code, format, database), timeout=60)
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

        return self._store(content, file, _is_binary)

    def _store(
        self, content: str | bytes, file: Path | None, _is_binary: bool
    ) -> Path | io.BytesIO | io.StringIO:
        """Write downloaded content to the cache file, or wrap it in memory."""
        if file:
            mode = "wb+" if _is_binary else "w+"
            # write then rename, so another process reading the cache (such as a
            # parallel test worker) never sees a partially written file
            partial = file.with_name(f".{file.name}.{os.getpid()}.partial")
            with open(partial, mode) as f:
                f.write(content)
            os.replace(partial, file)
            return Path(file)
        else:
            if _is_binary:
                if not isinstance(content, bytes):
                    raise ValueError(
                        "Binary content is not bytes, please check your format."
                    )
                return io.BytesIO(content)
            else:
                if not isinstance(content, str):
                    raise ValueError(
                        "Text content is not str, please check your format."
                    )
                return io.StringIO(content)

    def _download_pubchem(self, code: str) -> Path | io.StringIO:
        """Download the 3D conformer of a PubChem compound as SDF.

        Parameters
        ----------
        code : str
            A PubChem compound CID (all digits, e.g. '2244') or a compound name
            (e.g. 'aspirin').

        Returns
        -------
        Path or io.StringIO
            The cached `.sdf` file, or the SDF text when caching is disabled.

        Raises
        ------
        FileDownloadPDBError
            If the compound is not found, or PubChem has no 3D conformer for it.
        """
        namespace = "cid" if code.isdigit() else "name"
        file = self.cache / f"{code}.sdf" if self.cache else None
        if file and file.exists():
            return file

        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
            f"{namespace}/{quote(code, safe='')}/record/SDF?record_type=3d"
        )
        r = requests.get(url, timeout=60)
        if r.status_code == 404:
            if "No CID found" in r.text:
                raise FileDownloadPDBError(
                    f"No PubChem compound found with the name '{code}'."
                )
            # PubChem returns the same 404 for a CID with no 3D record and for a
            # CID that does not exist
            raise FileDownloadPDBError(
                f"PubChem has no 3D conformer for '{code}', or no such compound. "
                "PubChem only computes 3D conformers for single molecules (not salts "
                "or mixtures) of at most 50 heavy atoms and 15 rotatable bonds, made "
                "only of H, C, N, O, F, Si, P, S, Cl, Br and I (so no metals)."
            )
        try:
            r.raise_for_status()
        except requests.HTTPError as e:
            raise FileDownloadPDBError(f"Error downloading from PubChem: {e}")

        # text content is stored as a file or wrapped in StringIO, never BytesIO
        return cast(Path | io.StringIO, self._store(r.text, file, _is_binary=False))

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
        else:
            raise ValueError(f"Database {database} not currently supported.")
