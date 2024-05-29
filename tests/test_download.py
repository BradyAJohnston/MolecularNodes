from .constants import codes, data_dir
import tempfile
from biotite.structure.io import load_structure
import biotite.database.rcsb as rcsb
from molecularnodes.io.retrieve import download, FileDownloadPDBError
import os
import io
import pytest
import molecularnodes as mn


# currently can't figure out downloading from other services
databases = ["rcsb"]


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
        mn.io.download("4ozs", format=format, cache=tempfile.TemporaryDirectory().name)
    )
    struc_biotite = load_structure(
        rcsb.fetch(
            "4ozs", format=format, target_path=tempfile.TemporaryDirectory().name
        )
    )
    assert struc_download == struc_biotite


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("database", databases)
@pytest.mark.parametrize("format", ["pdb", "cif"])
def test_fetch_with_cache(tmpdir, code, format, database):
    cache_dir = tmpdir.mkdir("cache")
    file = mn.io.download(code, format, cache=str(cache_dir), database=database)

    assert isinstance(file, str)
    assert os.path.isfile(file)
    assert file.endswith(f"{code}.{format}")

    with open(file, "r") as f:
        content = f.read()
    assert content.startswith(_filestart(format))


databases = ["rcsb"]  # currently can't figure out downloading from the pdbe


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("database", databases)
@pytest.mark.parametrize("format", ["pdb", "cif"])
def test_fetch_without_cache(tmpdir, code, format, database):
    file = mn.io.download(code, format, cache=None, database=database)

    assert isinstance(file, io.StringIO)
    content = file.getvalue()
    assert content.startswith(_filestart(format))


@pytest.mark.parametrize("database", databases)
def test_fetch_with_invalid_format(database):
    code = "4OZS"
    format = "xyz"

    with pytest.raises(ValueError):
        mn.io.download(code, format, cache=None, database=database)


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("database", databases)
@pytest.mark.parametrize("format", ["bcif"])
def test_fetch_with_binary_format(tmpdir, code, database, format):
    cache_dir = tmpdir.mkdir("cache")
    file = mn.io.download(code, format, cache=str(cache_dir), database=database)

    assert isinstance(file, str)
    assert os.path.isfile(file)
    assert file.endswith(f"{code}.{format}")

    if format == "bcif":
        start = b"\x83\xa7"

    with open(file, "rb") as f:
        content = f.read()
    assert content.startswith(start)


# , 'bcif')) # TODO bcif tests once supported
@pytest.mark.parametrize("format", ("cif", "pdb"))
@pytest.mark.parametrize("code", ("A0A5E8G9H8", "A0A5E8G9T8", "K4PA18"))
def test_alphafold_download(format: str, code: str) -> None:
    file = mn.io.download(
        code=code, format=format, database="alphafold", cache=data_dir
    )
    mol = mn.io.load(file)

    assert mol.array
