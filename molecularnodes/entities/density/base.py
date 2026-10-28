import os
from abc import ABCMeta
from pathlib import Path
from typing import Union
import databpy
import numpy as np
from nodebpy import geometry as g
from nodebpy.builder import GeometrySocket
from ... import blender as bl
from ...nodes.geometry import SetGeometryVolume
from ..base import EntityType, MolecularEntity, MolecularTree
from .annotations import DensityAnnotationManager


def _is_set_geometry_volume(node) -> bool:
    tree = getattr(node, "node_tree", None)
    return tree is not None and tree.name.startswith(SetGeometryVolume._name)


class DensityTree(MolecularTree):
    """
    Node tree of a density, which keeps the source volume with its styles.

    Branches end in the join as usual, but the join feeds the output through a
    `Set Geometry Volume` node, which stores the tree's input volume in the
    `MN/Volume` bundle of the output. Density styles turn the volume into a mesh;
    the bundle keeps the grid available to `Color Sample Grid` without rendering it.
    """

    @property
    def join(self) -> GeometrySocket:
        output = self.geometry
        links = output.socket.links
        from_socket = links[0].from_socket if links else None
        keep = from_socket.node if from_socket is not None else None
        if keep is not None and _is_set_geometry_volume(keep):
            upstream = keep.inputs[0].links
            upstream_socket = upstream[0].from_socket if upstream else None
            upstream_node = (
                upstream_socket.node if upstream_socket is not None else None
            )
            if (
                upstream_node is not None
                and upstream_node.bl_idname == "GeometryNodeJoinGeometry"
            ):
                return self._wrap(upstream_node.inputs[0])
            join = g.JoinGeometry()
            if upstream_socket is not None:
                self._wrap(upstream_socket) >> join
            join >> self._wrap(keep.inputs[0])
            return join.i.geometry

        join_input = super().join
        join_node = join_input.socket.node
        assert join_node is not None
        join_output = self._wrap(join_node.outputs[0])
        SetGeometryVolume(geometry=join_output, volume=self.atoms) >> output
        return join_input


class Density(MolecularEntity, metaclass=ABCMeta):
    """
    Abstract base class for molecular density objects.

    """

    _tree_class = DensityTree
    _entity_type = EntityType.DENSITY

    def __init__(self, file_path: Union[str, Path]):
        super().__init__()
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
