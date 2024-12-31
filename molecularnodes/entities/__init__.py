from . import molecule, trajectory
from .base import MolecularEntity, EntityType
from .density import MN_OT_Import_Map
from .trajectory.dna import MN_OT_Import_OxDNA_Trajectory
from .ensemble import CellPack, StarFile, Ensemble
from .ensemble.ui import MN_OT_Import_Cell_Pack, MN_OT_Import_Star_File
from .molecule import Molecule, PDB, BCIF, CIF, SDF
from .molecule.ui import fetch, load_local, parse
from .trajectory import Trajectory, OXDNA

CLASSES = (
    [
        MN_OT_Import_Cell_Pack,
        MN_OT_Import_Map,
        MN_OT_Import_OxDNA_Trajectory,
        MN_OT_Import_Star_File,
    ]
    + trajectory.CLASSES
    + molecule.CLASSES
)
