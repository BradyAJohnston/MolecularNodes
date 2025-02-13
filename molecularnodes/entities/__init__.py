from . import molecule, trajectory, interaction
from .density import MN_OT_Import_Map
from .trajectory.dna import MN_OT_Import_OxDNA_Trajectory
from .ensemble import CellPack
from .ensemble import StarFile
from .ensemble.ui import MN_OT_Import_Cell_Pack, MN_OT_Import_Star_File
from .molecule.pdb import PDB
from .molecule.pdbx import BCIF, CIF
from .molecule.sdf import SDF
from .molecule.ui import fetch, load_local, parse
from .trajectory.trajectory import Trajectory

CLASSES = (
    [
        MN_OT_Import_Cell_Pack,
        MN_OT_Import_Map,
        MN_OT_Import_OxDNA_Trajectory,
        MN_OT_Import_Star_File,
    ]
    + trajectory.CLASSES
    + molecule.CLASSES
    + interaction.CLASSES
)
