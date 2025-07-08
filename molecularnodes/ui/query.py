import json
from dataclasses import dataclass
from typing import Dict, List, Optional
import bpy
import requests


@dataclass
class PDBStructureInfo:
    """Dataclass containing useful information about a PDB structure"""

    entry_id: str
    title: Optional[str] = None
    experimental_method: Optional[str] = None
    resolution: Optional[float] = None
    deposition_date: Optional[str] = None
    release_date: Optional[str] = None
    authors: Optional[List[str]] = None
    organism_names: Optional[List[str]] = None
    structure_classification: Optional[str] = None
    molecular_weight: Optional[float] = None
    polymer_entity_count: Optional[int] = None
    chain_count: Optional[int] = None
    residue_count: Optional[int] = None
    atom_count: Optional[int] = None
    space_group: Optional[str] = None
    unit_cell: Optional[Dict[str, float]] = None
    r_work: Optional[float] = None
    r_free: Optional[float] = None
    doi: Optional[str] = None
    pubmed_id: Optional[str] = None
    deposited_modeled_polymer_monomer_count: Optional[int] = None

    def __str__(self) -> str:
        """String representation of the structure info"""
        lines = [f"PDB Entry: {self.entry_id}"]
        if self.title:
            lines.append(f"Title: {self.title}")
        if self.experimental_method:
            lines.append(f"Method: {self.experimental_method}")
        if self.resolution:
            lines.append(f"Resolution: {self.resolution} Å")
        if self.organism_names:
            lines.append(f"Organisms: {', '.join(self.organism_names)}")
        return "\n".join(lines)

    def as_layout(self, layout: bpy.types.UILayout):
        layout = layout.column(align=True)
        layout.label(text=f"PDB: {self.entry_id.upper()}", icon="FILE_TEXT")
        layout.label(text=self.title)
        if self.resolution:
            layout.label(text=f"Resolution: {self.resolution} Å")
        if self.experimental_method:
            layout.label(text=f"Method: {self.experimental_method}")

    def to_json(self) -> str:
        """Convert the structure info to a JSON string"""
        return json.dumps(self.__dict__)

    @classmethod
    def from_json(cls, json_str: str) -> "PDBStructureInfo":
        """Create a PDBStructureInfo instance from a JSON string"""
        data = json.loads(json_str)
        return cls(**data)


class PDBQueryError(Exception):
    """Custom exception for PDB query errors"""

    def __init__(self, message: str, status_code: Optional[int] = None):
        self.message = message
        self.status_code = status_code
        super().__init__(self.message)


def query_pdb_structure(entry_id: str, timeout: int = 30) -> PDBStructureInfo:
    """
    Query the RCSB PDB database for basic structure information.

    Args:
        entry_id: The PDB entry ID (e.g., '4HHB', '1ABC')
        timeout: Request timeout in seconds (default: 30)

    Returns:
        PDBStructureInfo: Dataclass containing structure information

    Raises:
        PDBQueryError: If the query fails or structure doesn't exist
        ValueError: If entry_id is invalid format
    """
    # Validate entry_id format
    if not entry_id or not isinstance(entry_id, str):
        raise ValueError("entry_id must be a non-empty string")

    entry_id = entry_id.upper().strip()
    if len(entry_id) != 4:
        raise ValueError("entry_id must be exactly 4 characters long")

    # Construct the API URL
    url = f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}"

    try:
        # Make the API request
        response = requests.get(url, timeout=timeout)

        # Handle HTTP errors
        if response.status_code == 404:
            raise PDBQueryError(f"PDB entry '{entry_id}' not found", 404)
        elif response.status_code != 200:
            raise PDBQueryError(
                f"API request failed with status {response.status_code}: {response.text}",
                response.status_code,
            )

        # Parse JSON response
        try:
            data = response.json()
        except json.JSONDecodeError as e:
            raise PDBQueryError(f"Failed to parse JSON response: {e}")

        # Extract information from the response
        structure_info = PDBStructureInfo(entry_id=entry_id)

        # Basic structure information
        if "struct" in data:
            structure_info.title = data["struct"].get("title")

        # Experimental method
        if "exptl" in data and data["exptl"]:
            methods = [exp.get("method") for exp in data["exptl"] if exp.get("method")]
            structure_info.experimental_method = ", ".join(methods) if methods else None

        # Resolution
        if "rcsb_entry_info" in data:
            entry_info = data["rcsb_entry_info"]
            structure_info.resolution = entry_info.get("resolution_combined", [None])[0]
            structure_info.polymer_entity_count = (
                entry_info.get("polymer_entity_count_protein", 0)
                + entry_info.get("polymer_entity_count_DNA", 0)
                + entry_info.get("polymer_entity_count_RNA", 0)
            )
            structure_info.molecular_weight = entry_info.get("molecular_weight")
            structure_info.deposited_modeled_polymer_monomer_count = entry_info.get(
                "deposited_modeled_polymer_monomer_count"
            )

        # Dates
        if "pdbx_database_status" in data:
            db_status = data["pdbx_database_status"]
            structure_info.deposition_date = db_status.get(
                "recvd_initial_deposition_date"
            )
            structure_info.release_date = db_status.get("status_code_mr")

        # Authors
        if "audit_author" in data:
            authors = [
                author.get("name")
                for author in data["audit_author"]
                if author.get("name")
            ]
            structure_info.authors = authors if authors else None

        # Organism information
        if "rcsb_entity_source_organism" in data:
            organisms = []
            for org_data in data["rcsb_entity_source_organism"]:
                if "ncbi_scientific_name" in org_data:
                    organisms.append(org_data["ncbi_scientific_name"])
            structure_info.organism_names = organisms if organisms else None

        # Structure classification
        if "struct_keywords" in data:
            structure_info.structure_classification = data["struct_keywords"].get(
                "pdbx_keywords"
            )

        # Crystallographic information
        if "cell" in data:
            cell_data = data["cell"]
            structure_info.unit_cell = {
                "a": cell_data.get("length_a"),
                "b": cell_data.get("length_b"),
                "c": cell_data.get("length_c"),
                "alpha": cell_data.get("angle_alpha"),
                "beta": cell_data.get("angle_beta"),
                "gamma": cell_data.get("angle_gamma"),
            }

        if "symmetry" in data:
            structure_info.space_group = data["symmetry"].get("space_group_name_H-M")

        # Refinement statistics
        if "refine" in data and data["refine"]:
            refine_data = data["refine"][0]  # Usually one refinement entry
            structure_info.r_work = refine_data.get("ls_R_factor_R_work")
            structure_info.r_free = refine_data.get("ls_R_factor_R_free")

        # Citation information
        if "citation" in data:
            for citation in data["citation"]:
                if citation.get("id") == "primary":
                    structure_info.doi = citation.get("pdbx_database_id_DOI")
                    structure_info.pubmed_id = citation.get("pdbx_database_id_PubMed")
                    break

        # Additional counts
        if "rcsb_entry_info" in data:
            entry_info = data["rcsb_entry_info"]
            structure_info.chain_count = entry_info.get(
                "deposited_polymer_entity_instance_count"
            )
            structure_info.residue_count = entry_info.get(
                "deposited_modeled_polymer_monomer_count"
            )
            structure_info.atom_count = entry_info.get("deposited_atom_count")

        return structure_info

    except requests.exceptions.Timeout:
        raise PDBQueryError(f"Request timed out after {timeout} seconds")
    except requests.exceptions.ConnectionError:
        raise PDBQueryError("Failed to connect to RCSB PDB API")
    except requests.exceptions.RequestException as e:
        raise PDBQueryError(f"Request failed: {e}")


def _update_structure_display_info(self, context):
    code = context.scene.mn.import_code_pdb
    try:
        info = query_pdb_structure(code)
        context.scene.mn.import_display_info = info.to_json()
    except (PDBQueryError, ValueError) as e:
        context.scene.mn.import_display_info = f"ERROR: {e}"
