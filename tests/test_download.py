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
from molecularnodes.download import (
    CACHE_DIR,
    FileDownloadPDBError,
    StructureDownloader,
)
from .constants import codes

# currently can't figure out downloading from other services
DATABASES = ["rcsb"]
SLEEP_TIME = 2.0


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
    time.sleep(SLEEP_TIME)


def test_fail_download_pdb_large_structure_raises():
    downloader = StructureDownloader()
    with pytest.raises(FileDownloadPDBError):
        downloader.download("7D6Z", format="pdb")
    time.sleep(SLEEP_TIME)


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
            time.sleep(SLEEP_TIME)


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
    time.sleep(SLEEP_TIME)


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
    time.sleep(SLEEP_TIME)


DATABASES = ["rcsb"]  # currently can't figure out downloading from the pdbe


@pytest.mark.parametrize("database", DATABASES)
def test_fetch_with_invalid_format(database):
    code = "4OZS"
    format = "xyz"
    downloader = StructureDownloader(cache=None)

    with pytest.raises(ValueError):
        downloader.download(code, format, database=database)
    time.sleep(SLEEP_TIME)


# tests are failling intermittently - just ignoring for now TODO: stop ignoring
# limit to just 1 download test TODO: use API key to up requests in tests
# @pytest.mark.parametrize("format", ("cif",))
# # @pytest.mark.parametrize("format", ("pdb", "cif", "bcif"))
# @pytest.mark.parametrize(
#     "code",
#     (
#         "K4PA18",
#         # "G1JSI4",
#     ),
# )
# def test_fetch_alphafold(format: str, code: str, tmpdir) -> None:
#     time.sleep(SLEEP_TIME)
#     mol = mn.Molecule.fetch(code, format=format, database="alphafold", cache=tmpdir)
#     assert mol.array


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

    def test_url_alphafold_database_not_supported(self):
        # AlphaFold is now handled directly in download() using biotite.database.afdb
        # and should not be supported in _url()
        downloader = StructureDownloader(cache=None)
        with pytest.raises(
            ValueError, match="Database alphafold not currently supported"
        ):
            downloader._url("P12345", "cif", "alphafold")


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
