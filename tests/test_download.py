import gzip
import io
import os
import tempfile
import time
from pathlib import Path
from unittest.mock import Mock, patch
import biotite.database.rcsb as rcsb
import pytest
import requests
from biotite.structure.io import load_structure
import molecularnodes as mn
from molecularnodes.download import (
    CACHE_DIR,
    FileDownloadPDBError,
    StructureDownloader,
    get_alphafold_url,
)
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
    with tempfile.TemporaryDirectory() as temp_dir:
        downloader = StructureDownloader(cache=temp_dir)
        struc_download = load_structure(downloader.download("4ozs", format=format))

        with tempfile.TemporaryDirectory() as biotite_temp_dir:
            struc_biotite = load_structure(
                rcsb.fetch("4ozs", format=format, target_path=biotite_temp_dir)
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


# Test StructureDownloader initialization
class TestStructureDownloaderInit:
    def test_init_with_default_cache(self):
        downloader = StructureDownloader()
        assert downloader.cache == CACHE_DIR
        assert downloader.cache.exists()

    def test_init_with_custom_cache_string(self, tmpdir):
        cache_path = str(tmpdir.mkdir("custom_cache"))
        downloader = StructureDownloader(cache=cache_path)
        assert downloader.cache == Path(cache_path).absolute()
        assert downloader.cache.exists()

    def test_init_with_custom_cache_path(self, tmpdir):
        cache_path = tmpdir.mkdir("custom_cache")
        downloader = StructureDownloader(cache=Path(cache_path))
        assert downloader.cache == Path(cache_path).absolute()

    def test_init_with_no_cache(self):
        downloader = StructureDownloader(cache=None)
        assert downloader.cache is None

    def test_init_creates_cache_directory(self, tmpdir):
        nonexistent_path = tmpdir.join("new_cache")
        assert not nonexistent_path.exists()
        downloader = StructureDownloader(cache=str(nonexistent_path))
        assert downloader.cache.exists()


# Test URL generation
class TestUrlGeneration:
    def test_url_rcsb_cif(self):
        downloader = StructureDownloader(cache=None)
        url = downloader._url("1abc", "cif", "rcsb")
        assert url == "https://files.rcsb.org/download/1abc.cif"

    def test_url_rcsb_pdb(self):
        downloader = StructureDownloader(cache=None)
        url = downloader._url("1abc", "pdb", "rcsb")
        assert url == "https://files.rcsb.org/download/1abc.pdb"

    def test_url_rcsb_bcif(self):
        downloader = StructureDownloader(cache=None)
        url = downloader._url("1abc", "bcif", "rcsb")
        assert url == "https://models.rcsb.org/1abc.bcif"

    def test_url_pdb_database(self):
        downloader = StructureDownloader(cache=None)
        url = downloader._url("1abc", "cif", "pdb")
        assert url == "https://files.rcsb.org/download/1abc.cif"

    def test_url_wwpdb_database(self):
        downloader = StructureDownloader(cache=None)
        url = downloader._url("1abc", "cif", "wwpdb")
        assert url == "https://files.rcsb.org/download/1abc.cif"

    def test_url_unsupported_database(self):
        downloader = StructureDownloader(cache=None)
        with pytest.raises(
            ValueError, match="Database unsupported not currently supported"
        ):
            downloader._url("1abc", "cif", "unsupported")

    @patch("molecularnodes.download.get_alphafold_url")
    def test_url_alphafold_database(self, mock_get_alphafold_url):
        mock_get_alphafold_url.return_value = "https://alphafold.example.com/test.cif"
        downloader = StructureDownloader(cache=None)
        url = downloader._url("P12345", "cif", "alphafold")
        mock_get_alphafold_url.assert_called_once_with("P12345", "cif")
        assert url == "https://alphafold.example.com/test.cif"


# Test AlphaFold URL generation
class TestAlphaFoldUrl:
    @patch("requests.get")
    def test_get_alphafold_url_pdb(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = [
            {"pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.pdb"}
        ]
        mock_get.return_value = mock_response

        url = get_alphafold_url("P12345", "pdb")
        assert url == "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.pdb"
        mock_get.assert_called_once_with(
            "https://alphafold.ebi.ac.uk/api/prediction/P12345"
        )

    @patch("requests.get")
    def test_get_alphafold_url_cif(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = [
            {"cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.cif"}
        ]
        mock_get.return_value = mock_response

        url = get_alphafold_url("P12345", "cif")
        assert url == "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.cif"

    def test_get_alphafold_url_unsupported_format(self):
        with pytest.raises(
            ValueError,
            match="Format xyz not currently supported from AlphaFold database",
        ):
            get_alphafold_url("P12345", "xyz")

    @patch("requests.get")
    def test_get_alphafold_url_http_error(self, mock_get):
        mock_get.side_effect = requests.HTTPError("404 Not Found")
        with pytest.raises(requests.HTTPError):
            get_alphafold_url("INVALID", "pdb")


# Test cache behavior
class TestCacheBehavior:
    def test_cache_hit_returns_existing_file(self, tmpdir):
        cache_dir = tmpdir.mkdir("cache")
        existing_file = cache_dir.join("test.cif")
        existing_file.write("existing content")

        downloader = StructureDownloader(cache=str(cache_dir))

        with patch.object(downloader, "_url"), patch("requests.get") as mock_get:
            result = downloader.download("test", "cif")

            # Should return existing file without making HTTP request
            assert isinstance(result, Path)
            assert result.name == "test.cif"
            assert result.read_text() == "existing content"
            mock_get.assert_not_called()

    @patch("requests.get")
    def test_no_cache_returns_string_io(self, mock_get):
        mock_response = Mock()
        mock_response.text = "test content"
        mock_response.content = b"test content"
        mock_get.return_value = mock_response

        downloader = StructureDownloader(cache=None)
        result = downloader.download("test", "cif")

        assert isinstance(result, io.StringIO)
        assert result.getvalue() == "test content"

    @patch("requests.get")
    def test_no_cache_returns_bytes_io_for_binary(self, mock_get):
        mock_response = Mock()
        mock_response.content = b"binary content"
        mock_get.return_value = mock_response

        downloader = StructureDownloader(cache=None)
        result = downloader.download("test", "bcif")

        assert isinstance(result, io.BytesIO)
        assert result.getvalue() == b"binary content"


# Test error handling
class TestErrorHandling:
    @patch("requests.get")
    def test_http_error_raises_file_download_error(self, mock_get):
        mock_get.side_effect = requests.HTTPError("404 Not Found")

        downloader = StructureDownloader(cache=None)
        with pytest.raises(FileDownloadPDBError, match="404 Not Found"):
            downloader.download("nonexistent", "cif")

    def test_new_format_pdb_code_with_pdb_format_raises_error(self):
        downloader = StructureDownloader(cache=None)
        with pytest.raises(
            ValueError,
            match="New format PDB codes .* are not compatible with .pdb format",
        ):
            downloader.download("pdb_00001abc", "pdb")

    @patch("requests.get")
    def test_invalid_content_type_for_binary(self, mock_get):
        mock_response = Mock()
        mock_response.content = "not bytes"  # Should be bytes for binary format
        mock_get.return_value = mock_response

        downloader = StructureDownloader(cache=None)
        with pytest.raises(ValueError, match="Binary content is not bytes"):
            downloader.download("test", "bcif")

    @patch("requests.get")
    def test_invalid_content_type_for_text(self, mock_get):
        mock_response = Mock()
        mock_response.text = b"not string"  # Should be string for text format
        mock_get.return_value = mock_response

        downloader = StructureDownloader(cache=None)
        with pytest.raises(ValueError, match="Text content is not str"):
            downloader.download("test", "cif")


# Test gzip decompression
class TestGzipDecompression:
    @patch("requests.get")
    def test_gzipped_content_is_decompressed(self, mock_get):
        original_content = b"test binary content"
        gzipped_content = gzip.compress(original_content)

        mock_response = Mock()
        mock_response.content = gzipped_content
        mock_get.return_value = mock_response

        downloader = StructureDownloader(cache=None)
        result = downloader.download("test", "bcif")

        assert isinstance(result, io.BytesIO)
        assert result.getvalue() == original_content

    @patch("requests.get")
    def test_non_gzipped_binary_content_unchanged(self, mock_get):
        content = b"test binary content"

        mock_response = Mock()
        mock_response.content = content
        mock_get.return_value = mock_response

        downloader = StructureDownloader(cache=None)
        result = downloader.download("test", "bcif")

        assert isinstance(result, io.BytesIO)
        assert result.getvalue() == content


# Test input sanitization
class TestInputSanitization:
    @patch("requests.get")
    def test_code_with_whitespace_is_stripped(self, mock_get):
        mock_response = Mock()
        mock_response.text = "content"
        mock_get.return_value = mock_response

        downloader = StructureDownloader(cache=None)

        with patch.object(downloader, "_url") as mock_url:
            mock_url.return_value = "http://example.com"
            downloader.download("  test  ", "cif")
            # Check that the stripped code was used
            mock_url.assert_called_with("test", "cif", "rcsb")

    @patch("requests.get")
    def test_format_with_dot_is_stripped(self, mock_get):
        mock_response = Mock()
        mock_response.text = "content"
        mock_get.return_value = mock_response

        downloader = StructureDownloader(cache=None)

        with patch.object(downloader, "_url") as mock_url:
            mock_url.return_value = "http://example.com"
            downloader.download("test", ".cif")
            # Check that the stripped format was used
            mock_url.assert_called_with("test", "cif", "rcsb")


# Test file writing to cache
class TestFileWriting:
    @patch("requests.get")
    def test_text_file_written_to_cache(self, mock_get, tmpdir):
        content = "test content"
        mock_response = Mock()
        mock_response.text = content
        mock_get.return_value = mock_response

        cache_dir = tmpdir.mkdir("cache")
        downloader = StructureDownloader(cache=str(cache_dir))

        result = downloader.download("test", "cif")

        assert isinstance(result, Path)
        assert result.read_text() == content

    @patch("requests.get")
    def test_binary_file_written_to_cache(self, mock_get, tmpdir):
        content = b"binary content"
        mock_response = Mock()
        mock_response.content = content
        mock_get.return_value = mock_response

        cache_dir = tmpdir.mkdir("cache")
        downloader = StructureDownloader(cache=str(cache_dir))

        result = downloader.download("test", "bcif")

        assert isinstance(result, Path)
        assert result.read_bytes() == content
