import os
from abc import ABCMeta
from pathlib import Path
from typing import Union
import databpy
import numpy as np
from ... import blender as bl
from ..base import EntityType, MolecularEntity


class Density(MolecularEntity, metaclass=ABCMeta):
    """
    Abstract base class for molecular density objects.

    """

    def __init__(self, file_path: Union[str, Path]):
        super().__init__()
        self._entity_type = EntityType.DENSITY
        self.file_path: Path = bl.path_resolve(file_path)
        self.grid = None
        self.file_vdb: str
        self.threshold: float

    def named_attribute(self, name: str, evaluate: bool = True) -> np.ndarray:
        obj = bl.mesh.evaluate_using_mesh(self.object)
        return databpy.named_attribute(obj, name, evaluate=True)

    def path_to_vdb(self, file: str, center: bool = False, invert: bool = False):
        """
        Convert a file path to a corresponding VDB file path.

        Parameters
        ----------
        file : str
            The path of the original file.
        center : bool, default=False
            If True, the density will be centered at the origin.
        invert : bool, default=False
            If True, the density values will be inverted.

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
