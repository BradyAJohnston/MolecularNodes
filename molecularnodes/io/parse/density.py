from abc import ABCMeta
import os
import bpy

from typing import Union, Optional
from pathlib import Path


class Density(metaclass=ABCMeta):
    """
    Abstract base class for molecular density objects.

    """

    def __init__(self, file_path: Union[str, Path]) -> None:
        self.file_path = file_path
        self.grid = None
        self.file_vdb: Union[Path, str]
        self.threshold: float = 1.0
        self.object: Optional[bpy.types.Object] = None

    def path_to_vdb(self, file: Union[Path, str], center: False, invert: False) -> Path:
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
        return Path(file_path)
