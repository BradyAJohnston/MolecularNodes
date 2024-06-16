from .cellpack import MN_OT_Import_Cell_Pack
from .density import MN_OT_Import_Map
from .dna import MN_OT_Import_OxDNA_Trajectory
from .local import MN_OT_Import_Protein_Local, load
from . import md
from .alphafold import MN_OT_Import_AlphaFold
from .parse import BCIF, CIF, PDB, SDF, CellPack, StarFile, MNUniverse
from .retrieve import download
from .star import MN_OT_Import_Star_File
from .wwpdb import MN_OT_Import_wwPDB, fetch

ops_io = [
    MN_OT_Import_AlphaFold,
    MN_OT_Import_Cell_Pack,
    MN_OT_Import_Map,
    MN_OT_Import_OxDNA_Trajectory,
    MN_OT_Import_Protein_Local,
    MN_OT_Import_Star_File,
    MN_OT_Import_wwPDB,
] + md.CLASSES
