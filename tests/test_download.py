import io
import os
import tempfile

import biotite.database.rcsb as rcsb
import pytest
from biotite.structure.io import load_structure

import molecularnodes as mn
from molecularnodes.download import FileDownloadPDBError, download

from .constants import codes

# currently can't figure out downloading from other services
DATABASES = ["rcsb"]


def _filestart(format):
    if format == "cif":
        return "data_"
    else:
        return "HEADER"


def test_download_raises_error_on_invalid_format():
    with pytest.raises(ValueError) as excinfo:
        download("1abc", "invalid_format")
    assert (
        "File format 'invalid_format' not in: supported_formats=['cif', 'pdb', 'bcif']"
        in str(excinfo.value)
    )


def test_fail_download_pdb_large_structure_raises():
    with pytest.raises(FileDownloadPDBError) as excinfo:
        download("7D6Z", format="pdb")

    assert (
        "There was an error downloading the file from the Protein Data Bank. PDB or format for PDB code may not be available."
        in str(excinfo.value)
    )


@pytest.mark.parametrize("format", ["cif", "bcif", "pdb"])
def test_compare_biotite(format):
    struc_download = load_structure(
        download("4ozs", format=format, cache=tempfile.TemporaryDirectory().name)
    )
    struc_biotite = load_structure(
        rcsb.fetch(
            "4ozs", format=format, target_path=tempfile.TemporaryDirectory().name
        )
    )
    assert struc_download == struc_biotite


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("database", DATABASES)
@pytest.mark.parametrize("format", ["pdb", "cif"])
def test_fetch_with_cache(tmpdir, code, format, database):
    cache_dir = tmpdir.mkdir("cache")
    file = download(code, format, cache=str(cache_dir), database=database)

    assert isinstance(file, str)
    assert os.path.isfile(file)
    assert file.endswith(f"{code}.{format}")

    with open(file, "r") as f:
        content = f.read()
    assert content.startswith(_filestart(format))


DATABASES = ["rcsb"]  # currently can't figure out downloading from the pdbe


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("database", DATABASES)
@pytest.mark.parametrize("format", ["pdb", "cif"])
def test_fetch_without_cache(tmpdir, code, format, database):
    file = download(code, format, cache=None, database=database)

    assert isinstance(file, io.StringIO)
    content = file.getvalue()
    assert content.startswith(_filestart(format))


@pytest.mark.parametrize("database", DATABASES)
def test_fetch_with_invalid_format(database):
    code = "4OZS"
    format = "xyz"

    with pytest.raises(ValueError):
        download(code, format, cache=None, database=database)


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("database", DATABASES)
@pytest.mark.parametrize("format", ["bcif"])
def test_fetch_with_binary_format(tmpdir, code, database, format):
    cache_dir = tmpdir.mkdir("cache")
    file = download(code, format, cache=str(cache_dir), database=database)

    assert isinstance(file, str)
    assert os.path.isfile(file)
    assert file.endswith(f"{code}.{format}")

    if format == "bcif":
        start = b"\x83\xa7"

    with open(file, "rb") as f:
        content = f.read()
    assert content.startswith(start)


# TODO BCIF is supported elsewhere in the package but can't currently be parsed properly
# I think there is something weird going on with the alphafold formatted bcif files


@pytest.mark.parametrize("format", ("cif", "pdb"))
@pytest.mark.parametrize("code", ("A0A5E8G9H8", "A0A5E8G9T8", "K4PA18"))
def test_alphafold_download(format: str, code: str, tmpdir) -> None:
    file = download(code=code, format=format, database="alphafold", cache=tmpdir)

    mol = mn.entities.load_local(file)

    assert mol.array
