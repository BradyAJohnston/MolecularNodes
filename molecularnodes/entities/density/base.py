import os
from abc import ABCMeta
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Union
import databpy
import numpy as np
from nodebpy.builder import GeometrySocket
from ... import blender as bl
from ..base import EntityType, MolecularEntity, MolecularTree, ResetSockets
from .annotations import DensityAnnotationManager


class DensityTree(MolecularTree):
    """
    Node tree of a density.

    The density's object is an empty mesh, so the tree has no input geometry:
    each branch starts with an `Import Volume` node that reads the density's
    ``.vdb`` files, and ends in the join.

    >>> with density.tree.reset() as (_, join):
    ...     ImportVolume(path) >> DensityStyleSurface() >> join
    """

    @property
    def atoms(self) -> GeometrySocket:
        raise AttributeError(
            "Density trees have no input geometry; start each branch with an "
            "`ImportVolume` node"
        )

    @contextmanager
    def reset(self, output: str = "Geometry") -> Iterator[ResetSockets]:
        """
        Clear the tree and build within it.

        Yields ``(None, join)``: there is no input geometry to build from, so
        branches start with an `Import Volume` node and end in ``join``.
        """
        with self:
            self.clear()
            self.outputs.geometry(output)
            yield ResetSockets(None, self.join)


class Density(MolecularEntity, metaclass=ABCMeta):
    """
    Abstract base class for molecular density objects.

    """

    _tree_class = DensityTree

    def __init__(self, file_path: Union[str, Path]):
        super().__init__()
        self._entity_type = EntityType.DENSITY
        self.file_path: Path = bl.path_resolve(file_path)
        self.file_vdb: str
        self.annotations = DensityAnnotationManager(self)

    def named_attribute(self, name: str, evaluate: bool = True) -> np.ndarray:
        # volume objects can't be evaluated to a mesh data-block directly, but the
        # evaluated geometry holds the mesh their node tree generates
        return databpy.GeometrySet(self.object).named_attribute(name)

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

    def __getstate__(self):
        """Custom serialization."""
        state = super().__getstate__()
        # Remove objects with circular references or PyCapsules
        if "annotations" in state:
            del state["annotations"]
        return state

    def __setstate__(self, state):
        """Custom deserialization."""
        self.__dict__.update(state)
        # Recreate objects with circular references
        if not hasattr(self, "annotations"):
            self.annotations = DensityAnnotationManager(self)
