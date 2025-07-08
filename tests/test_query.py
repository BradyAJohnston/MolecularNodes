import json
from unittest.mock import Mock, patch
import pytest
from requests.exceptions import ConnectionError, RequestException, Timeout
from molecularnodes.ui.query import (
    PDBQueryError,
    PDBStructureInfo,
    query_pdb_structure,
)

# Import the functions and classes to test
# from rcsb_query import query_pdb_structure, PDBStructureInfo, PDBQueryError


class TestPDBStructureInfo:
    """Test the PDBStructureInfo dataclass"""

    def test_dataclass_initialization(self):
        """Test basic initialization of PDBStructureInfo"""
        info = PDBStructureInfo(entry_id="4HHB")
        assert info.entry_id == "4HHB"
        assert info.title is None
        assert info.experimental_method is None
        assert info.resolution is None

    def test_dataclass_with_data(self):
        """Test PDBStructureInfo with actual data"""
        info = PDBStructureInfo(
            entry_id="4HHB",
            title="Test Structure",
            experimental_method="X-RAY DIFFRACTION",
            resolution=2.5,
            authors=["Smith, J.", "Doe, J."],
            organism_names=["Homo sapiens"],
        )
        assert info.entry_id == "4HHB"
        assert info.title == "Test Structure"
        assert info.experimental_method == "X-RAY DIFFRACTION"
        assert info.resolution == 2.5
        assert len(info.authors) == 2
        assert info.organism_names[0] == "Homo sapiens"

    def test_string_representation(self):
        """Test __str__ method"""
        info = PDBStructureInfo(
            entry_id="4HHB",
            title="Test Structure",
            experimental_method="X-RAY DIFFRACTION",
            resolution=2.5,
            organism_names=["Homo sapiens"],
        )
        str_repr = str(info)
        assert "PDB Entry: 4HHB" in str_repr
        assert "Title: Test Structure" in str_repr
        assert "Method: X-RAY DIFFRACTION" in str_repr
        assert "Resolution: 2.5 Å" in str_repr
        assert "Organisms: Homo sapiens" in str_repr


class TestPDBQueryError:
    """Test the PDBQueryError exception"""

    def test_basic_error(self):
        """Test basic error creation"""
        error = PDBQueryError("Test error message")
        assert error.message == "Test error message"
        assert error.status_code is None
        assert str(error) == "Test error message"

    def test_error_with_status_code(self):
        """Test error with status code"""
        error = PDBQueryError("Not found", 404)
        assert error.message == "Not found"
        assert error.status_code == 404


class TestQueryPDBStructure:
    """Test the query_pdb_structure function"""

    def test_invalid_entry_id_formats(self):
        """Test validation of entry_id formats"""
        # Test empty string
        with pytest.raises(ValueError, match="entry_id must be a non-empty string"):
            query_pdb_structure("")

        # Test None
        with pytest.raises(ValueError, match="entry_id must be a non-empty string"):
            query_pdb_structure(None)

        # Test wrong length
        with pytest.raises(
            ValueError, match="entry_id must be exactly 4 characters long"
        ):
            query_pdb_structure("123")

        with pytest.raises(
            ValueError, match="entry_id must be exactly 4 characters long"
        ):
            query_pdb_structure("12345")

    def test_entry_id_normalization(self):
        """Test that entry_id is properly normalized (uppercase, stripped)"""
        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.json.return_value = {"struct": {"title": "Test"}}
            mock_get.return_value = mock_response

            # Test lowercase conversion
            query_pdb_structure("4hhb")
            mock_get.assert_called_with(
                "https://data.rcsb.org/rest/v1/core/entry/4HHB", timeout=30
            )

            # Test whitespace stripping
            query_pdb_structure(" 4HHB ")
            mock_get.assert_called_with(
                "https://data.rcsb.org/rest/v1/core/entry/4HHB", timeout=30
            )

    @patch("requests.get")
    def test_successful_query(self, mock_get):
        """Test successful API query with mock data"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "struct": {"title": "Test Hemoglobin Structure"},
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {
                "resolution_combined": [2.5],
                "polymer_entity_count_protein": 4,
                "polymer_entity_count_DNA": 0,
                "polymer_entity_count_RNA": 0,
                "molecular_weight": 64500.0,
                "deposited_modeled_polymer_monomer_count": 574,
            },
            "audit_author": [{"name": "Perutz, M.F."}, {"name": "Rossmann, M.G."}],
            "rcsb_entity_source_organism": [{"ncbi_scientific_name": "Homo sapiens"}],
            "struct_keywords": {"pdbx_keywords": "OXYGEN TRANSPORT"},
            "cell": {
                "length_a": 63.150,
                "length_b": 83.590,
                "length_c": 53.800,
                "angle_alpha": 90.00,
                "angle_beta": 99.34,
                "angle_gamma": 90.00,
            },
            "symmetry": {"space_group_name_H-M": "P 1 21 1"},
            "refine": [{"ls_R_factor_R_work": 0.160, "ls_R_factor_R_free": 0.230}],
        }
        mock_get.return_value = mock_response

        result = query_pdb_structure("4HHB")

        assert result.entry_id == "4HHB"
        assert result.title == "Test Hemoglobin Structure"
        assert result.experimental_method == "X-RAY DIFFRACTION"
        assert result.resolution == 2.5
        assert result.polymer_entity_count == 4
        assert result.authors == ["Perutz, M.F.", "Rossmann, M.G."]
        assert result.organism_names == ["Homo sapiens"]
        assert result.structure_classification == "OXYGEN TRANSPORT"
        assert result.space_group == "P 1 21 1"
        assert result.r_work == 0.160
        assert result.r_free == 0.230
        assert result.unit_cell["a"] == 63.150
        assert result.unit_cell["beta"] == 99.34

    @patch("requests.get")
    def test_404_not_found(self, mock_get):
        """Test handling of 404 response"""
        mock_response = Mock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response

        with pytest.raises(PDBQueryError) as exc_info:
            query_pdb_structure("XXXX")

        assert "PDB entry 'XXXX' not found" in str(exc_info.value)
        assert exc_info.value.status_code == 404

    @patch("requests.get")
    def test_other_http_errors(self, mock_get):
        """Test handling of other HTTP errors"""
        mock_response = Mock()
        mock_response.status_code = 500
        mock_response.text = "Internal Server Error"
        mock_get.return_value = mock_response

        with pytest.raises(PDBQueryError) as exc_info:
            query_pdb_structure("4HHB")

        assert "API request failed with status 500" in str(exc_info.value)
        assert exc_info.value.status_code == 500

    @patch("requests.get")
    def test_json_decode_error(self, mock_get):
        """Test handling of invalid JSON response"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.side_effect = json.JSONDecodeError("Invalid JSON", "", 0)
        mock_get.return_value = mock_response

        with pytest.raises(PDBQueryError) as exc_info:
            query_pdb_structure("4HHB")

        assert "Failed to parse JSON response" in str(exc_info.value)

    @patch("requests.get")
    def test_timeout_error(self, mock_get):
        """Test handling of timeout"""
        mock_get.side_effect = Timeout()

        with pytest.raises(PDBQueryError) as exc_info:
            query_pdb_structure("4HHB", timeout=5)

        assert "Request timed out after 5 seconds" in str(exc_info.value)

    @patch("requests.get")
    def test_connection_error(self, mock_get):
        """Test handling of connection error"""
        mock_get.side_effect = ConnectionError()

        with pytest.raises(PDBQueryError) as exc_info:
            query_pdb_structure("4HHB")

        assert "Failed to connect to RCSB PDB API" in str(exc_info.value)

    @patch("requests.get")
    def test_general_request_exception(self, mock_get):
        """Test handling of general request exceptions"""
        mock_get.side_effect = RequestException("Network error")

        with pytest.raises(PDBQueryError) as exc_info:
            query_pdb_structure("4HHB")

        assert "Request failed: Network error" in str(exc_info.value)

    @patch("requests.get")
    def test_minimal_response(self, mock_get):
        """Test handling of minimal API response"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {}  # Empty response
        mock_get.return_value = mock_response

        result = query_pdb_structure("4HHB")

        assert result.entry_id == "4HHB"
        assert result.title is None
        assert result.experimental_method is None
        assert result.resolution is None

    @patch("requests.get")
    def test_multiple_experimental_methods(self, mock_get):
        """Test handling of multiple experimental methods"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "exptl": [
                {"method": "X-RAY DIFFRACTION"},
                {"method": "NEUTRON DIFFRACTION"},
            ]
        }
        mock_get.return_value = mock_response

        result = query_pdb_structure("4HHB")

        assert result.experimental_method == "X-RAY DIFFRACTION, NEUTRON DIFFRACTION"

    @patch("requests.get")
    def test_custom_timeout(self, mock_get):
        """Test custom timeout parameter"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {}
        mock_get.return_value = mock_response

        query_pdb_structure("4HHB", timeout=60)

        mock_get.assert_called_with(
            "https://data.rcsb.org/rest/v1/core/entry/4HHB", timeout=60
        )

    @patch("requests.get")
    def test_citation_parsing(self, mock_get):
        """Test parsing of citation information"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "citation": [
                {
                    "id": "primary",
                    "pdbx_database_id_DOI": "10.1038/171740a0",
                    "pdbx_database_id_PubMed": "14882605",
                },
                {
                    "id": "secondary",
                    "pdbx_database_id_DOI": "10.1016/other",
                    "pdbx_database_id_PubMed": "12345678",
                },
            ]
        }
        mock_get.return_value = mock_response

        result = query_pdb_structure("4HHB")

        # Should only get the primary citation
        assert result.doi == "10.1038/171740a0"
        assert result.pubmed_id == "14882605"


class TestIntegration:
    """Integration tests (these will make actual API calls)"""

    @pytest.mark.integration
    def test_real_api_call_valid_structure(self):
        """Test with a real API call to a known structure"""
        # Using 1CRN which is a small, stable structure
        try:
            result = query_pdb_structure("1CRN")
            assert result.entry_id == "1CRN"
            assert result.title is not None
            assert result.experimental_method is not None
            # 1CRN should be NMR
            assert "NMR" in result.experimental_method.upper()
        except Exception as e:
            pytest.skip(f"Real API call failed: {e}")

    @pytest.mark.integration
    def test_real_api_call_invalid_structure(self):
        """Test with a real API call to an invalid structure"""
        with pytest.raises(PDBQueryError) as exc_info:
            query_pdb_structure("ZZZZ")

        assert exc_info.value.status_code == 404
        assert "not found" in str(exc_info.value).lower()


# Pytest fixtures
@pytest.fixture
def sample_pdb_data():
    """Sample PDB data for testing"""
    return {
        "struct": {"title": "Test Structure"},
        "exptl": [{"method": "X-RAY DIFFRACTION"}],
        "rcsb_entry_info": {
            "resolution_combined": [2.0],
            "polymer_entity_count_protein": 1,
            "polymer_entity_count_DNA": 0,
            "polymer_entity_count_RNA": 0,
        },
        "audit_author": [{"name": "Test Author"}],
        "rcsb_entity_source_organism": [{"ncbi_scientific_name": "Test organism"}],
    }


@pytest.fixture
def mock_successful_response(sample_pdb_data):
    """Mock successful API response"""
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.json.return_value = sample_pdb_data
    return mock_response


# Parametrized tests for edge cases
@pytest.mark.parametrize(
    "entry_id,expected_normalized",
    [("4hhb", "4HHB"), ("  4HHB  ", "4HHB"), ("1abc", "1ABC"), ("1A2B", "1A2B")],
)
def test_entry_id_normalization_parametrized(entry_id, expected_normalized):
    """Test entry ID normalization with various inputs"""
    with patch("requests.get") as mock_get:
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {}
        mock_get.return_value = mock_response

        query_pdb_structure(entry_id)
        expected_url = f"https://data.rcsb.org/rest/v1/core/entry/{expected_normalized}"
        mock_get.assert_called_with(expected_url, timeout=30)


if __name__ == "__main__":
    # Run tests with: python -m pytest test_rcsb_query.py -v
    # Run integration tests: python -m pytest test_rcsb_query.py -v -m integration
    pytest.main([__file__, "-v"])
