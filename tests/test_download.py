import os
import io
import pytest
import molecularnodes as mn
import biotite.database.rcsb as rcsb
from biotite.structure.io import load_structure
import tempfile

from .constants import codes

# currently can't figure out downloading from other services
databases = ['rcsb']


def _filestart(format):
    if format == "cif":
        return 'data_'
    else:
        return 'HEADER'


@pytest.mark.parametrize('format', ['cif', 'mmtf', 'pdb'])
def test_compare_biotite(format):
    struc_download = load_structure(mn.io.download(
        '4ozs', format=format, cache=tempfile.TemporaryDirectory().name))
    struc_biotite = load_structure(rcsb.fetch(
        '4ozs', format=format, target_path=tempfile.TemporaryDirectory().name))
    assert struc_download == struc_biotite


@pytest.mark.parametrize('code', codes)
@pytest.mark.parametrize('database', databases)
@pytest.mark.parametrize('format', ['pdb', 'cif'])
def test_fetch_with_cache(tmpdir, code, format, database):
    cache_dir = tmpdir.mkdir("cache")
    file = mn.io.download(code, format, cache=str(
        cache_dir), database=database)

    assert isinstance(file, str)
    assert os.path.isfile(file)
    assert file.endswith(f"{code}.{format}")

    with open(file, "r") as f:
        content = f.read()
    assert content.startswith(_filestart(format))


databases = ['rcsb']  # currently can't figure out downloading from the pdbe


@pytest.mark.parametrize('code', codes)
@pytest.mark.parametrize('database', databases)
@pytest.mark.parametrize('format', ['pdb', 'cif'])
def test_fetch_without_cache(tmpdir, code, format, database):
    file = mn.io.download(code, format, cache=None, database=database)

    assert isinstance(file, io.StringIO)
    content = file.getvalue()
    assert content.startswith(_filestart(format))


@pytest.mark.parametrize('database', databases)
def test_fetch_with_invalid_format(database):
    code = '4OZS'
    format = "xyz"

    with pytest.raises(ValueError):
        mn.io.download(code, format, cache=None, database=database)


@pytest.mark.parametrize('code', codes)
@pytest.mark.parametrize('database', databases)
@pytest.mark.parametrize('format', ['bcif', 'mmtf'])
def test_fetch_with_binary_format(tmpdir, code, database, format):
    cache_dir = tmpdir.mkdir("cache")
    file = mn.io.download(code, format, cache=str(
        cache_dir), database=database)

    assert isinstance(file, str)
    assert os.path.isfile(file)
    assert file.endswith(f"{code}.{format}")

    if format == "bcif":
        start = b"\x83\xa7"
    elif format == "mmtf":
        start = b"\xde\x00"

    with open(file, "rb") as f:
        content = f.read()
    assert content.startswith(start)


# , 'bcif')) # TODO bcif tests once supported
@pytest.mark.parametrize('format', ('cif', 'pdb'))
@pytest.mark.parametrize('code', ('A0A5E8G9H8', 'A0A5E8G9T8', 'K4PA18'))
def test_alphafold_download(format: str, code: str) -> None:
    file = mn.io.download(code=code, format=format, database='alphafold')
    if format == "bcif":
        model = mn.io.parse.BCIF(file)
    elif format == "cif":
        model = mn.io.parse.CIF(file)
    elif format == "pdb":
        model = mn.io.parse.PDB(file)

    assert model.array
