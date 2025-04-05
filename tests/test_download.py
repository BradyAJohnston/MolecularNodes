import io
import os
import tempfile
import time
from pathlib import Path
import biotite.database.rcsb as rcsb
import pytest
from biotite.structure.io import load_structure
import molecularnodes as mn
from molecularnodes.download import FileDownloadPDBError, StructureDownloader
from .constants import codes

# currently can't figure out downloading from other services
DATABASES = ["rcsb"]


def _filestart(format):
    if format == "cif":
        return "data_"
    else:
        return "HEADER"


def test_download_raises_error_on_invalid_format():
    downloader = StructureDownloader()
    with pytest.raises(ValueError) as excinfo:
        downloader.download("1abc", "invalid_format")
    assert (
        "File format 'invalid_format' not in: supported_formats=['cif', 'pdb', 'bcif']"
        in str(excinfo.value)
    )
    time.sleep(0.5)


def test_fail_download_pdb_large_structure_raises():
    downloader = StructureDownloader()
    with pytest.raises(FileDownloadPDBError):
        downloader.download("7D6Z", format="pdb")
    time.sleep(0.5)


@pytest.mark.parametrize("format", ["cif", "bcif", "pdb"])
def test_compare_biotite(format):
    downloader = StructureDownloader(cache=tempfile.TemporaryDirectory().name)
    struc_download = load_structure(downloader.download("4ozs", format=format))
    struc_biotite = load_structure(
        rcsb.fetch(
            "4ozs", format=format, target_path=tempfile.TemporaryDirectory().name
        )
    )
    assert struc_download == struc_biotite
    time.sleep(0.5)


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("database", DATABASES)
@pytest.mark.parametrize("format", ["pdb", "cif", "bcif"])
def test_fetch_with_cache(tmpdir, code, format, database):
    cache_dir = tmpdir.mkdir("cache")
    downloader = StructureDownloader(cache=str(cache_dir))
    file = downloader.download(code, format, database=database)

    assert isinstance(file, Path)
    assert os.path.isfile(file)
    assert file.name == f"{code}.{format}"

    if format == "bcif":
        # Binary format needs to be opened in binary mode
        with open(file, "rb") as f:
            content = f.read()
        # Use the binary file starts defined in test_fetch_with_binary_format
        assert content.startswith(b"\x83\xa7")
    else:
        # Text formats can be opened normally
        with open(file, "r") as f:
            content = f.read()
        assert content.startswith(_filestart(format))
    time.sleep(0.5)


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("database", DATABASES)
@pytest.mark.parametrize("format", ["pdb", "cif", "bcif"])
def test_fetch_new_code(tmpdir, code, format, database):
    cache_dir = tmpdir.mkdir("cache")
    downloader = StructureDownloader(cache=str(cache_dir))
    code = f"pdb_{code.upper().rjust(8, '0')}"
    if format == "pdb":
        with pytest.raises(ValueError):
            downloader.download(code, format, database=database)
        return
    if format == "bcif":
        with pytest.raises(FileDownloadPDBError):
            downloader.download(code, format, database=database)
        return
    file = downloader.download(code, format, database=database)

    assert isinstance(file, Path)
    assert os.path.isfile(file)
    assert file.name == f"{code}.{format}"

    with open(file, "r") as f:
        content = f.read()
    assert content.startswith(_filestart(format))
    time.sleep(0.5)


DATABASES = ["rcsb"]  # currently can't figure out downloading from the pdbe


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("database", DATABASES)
@pytest.mark.parametrize("format", ["pdb", "cif"])
def test_fetch_without_cache(tmpdir, code, format, database):
    downloader = StructureDownloader(cache=None)
    file = downloader.download(code, format, database=database)

    assert isinstance(file, io.StringIO)
    content = file.getvalue()
    assert content.startswith(_filestart(format))
    time.sleep(0.5)


@pytest.mark.parametrize("database", DATABASES)
def test_fetch_with_invalid_format(database):
    code = "4OZS"
    format = "xyz"
    downloader = StructureDownloader(cache=None)

    with pytest.raises(ValueError):
        downloader.download(code, format, database=database)
    time.sleep(0.5)


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("database", DATABASES)
@pytest.mark.parametrize("format", ["bcif"])
def test_fetch_with_binary_format(tmpdir, code, database, format):
    cache_dir = tmpdir.mkdir("cache")
    downloader = StructureDownloader(cache=str(cache_dir))
    file = downloader.download(code, format, database=database)

    assert isinstance(file, Path)
    assert os.path.isfile(file)
    assert file.name == f"{code}.{format}"

    start = {
        "bcif": b"\x83\xa7",
        "cif": b"data_",
        "pdb": b"HEADER",
    }[format]
    with open(file, "rb") as f:
        content = f.read()
    assert content.startswith(start)
    time.sleep(0.5)


# TODO BCIF is supported elsewhere in the package but can't currently be parsed properly
# I think there is something weird going on with the alphafold formatted bcif files


@pytest.mark.parametrize("format", ("cif", "pdb"))
@pytest.mark.parametrize("code", ("A0A5E8G9H8", "A0A5E8G9T8", "K4PA18"))
def test_alphafold_download(format: str, code: str, tmpdir) -> None:
    downloader = StructureDownloader(cache=tmpdir)
    file = downloader.download(code=code, format=format, database="alphafold")

    mol = mn.Molecule.load(file)

    assert mol.array
    time.sleep(0.5)
