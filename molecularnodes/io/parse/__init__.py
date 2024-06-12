"""
A subpackge which provides classes for parsing the different macromolecular data formats.
"""

from .pdbx import CIF, BCIF
# from .bcif import BCIF
# from .cif import CIF
from .pdb import PDB
from .cellpack import CellPack
from .star import StarFile, NDJSON
from .sdf import SDF
from .mda import MDAnalysisSession
from .mrc import MRC
