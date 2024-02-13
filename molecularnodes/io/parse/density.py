from abc import ABCMeta
import os
import numpy as np
import bpy


class Density(metaclass=ABCMeta):
    """
    Abstract base class for molecular density objects.

    """

    def __init__(self, file_path):
        self.file_path: str = None
        self.grid = None
        self.file_vdb: str = None
        self.threshold: float = None
        self.object: bpy.types.Object = None

    def path_to_vdb(self, file: str):
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
        file_name = name + '.vdb'
        file_path = os.path.join(folder_path, file_name)
        return file_path
