from abc import ABCMeta
import os
import numpy as np
import bpy
from typing import Union
from pathlib import Path
from ..entity import MolecularEntity, EntityType
from ... import blender as bl


class Density(MolecularEntity, metaclass=ABCMeta):
    """
    Abstract base class for molecular density objects.

    """

    def __init__(self, file_path: Union[str, Path]):
        self._entity_type = EntityType.DENSITY
        self.file_path: Path = bl.path_resolve(file_path)
        self.grid = None
        self.file_vdb: str
        self.threshold: float

    def path_to_vdb(self, file: str, center: False, invert: False):
        """
        Convert a file path to a corresponding VDB file path.

        Parameters
        ----------
        file : str
            The path of the original file.

        Returns
        -------
        str
            The path of the corresponding VDB file.
        """
        # Set up file paths
        folder_path = os.path.dirname(file)
        name = os.path.basename(file).split(".")[0]
        name += "_center" if center else ""
        name += "_invert" if invert else ""
        file_name = name + ".vdb"
        file_path = os.path.join(folder_path, file_name)
        return file_path
