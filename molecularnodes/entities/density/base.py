import os
from abc import ABCMeta
from contextlib import contextmanager
from pathlib import Path
from typing import Generator, Union, overload
import bpy
import databpy
import numpy as np
from nodebpy.builder import GeometrySocket
from ... import blender as bl
from ...nodes.geometry import ImportVolume
from ..base import EntityType, MolecularEntity, MolecularTree, ResetSockets
from .annotations import DensityAnnotationManager


class DensityTree(MolecularTree):
    """
    Node tree of a density.

    The density's object is an empty mesh, so the tree has no input geometry:
    each branch starts with an `Import Volume` node that reads the density's
    ``.vdb`` files, and ends in the join.

    >>> with density.tree.reset() as (volume, join):
    ...     volume >> DensityStyleSurface() >> join
    """

    @property
    def atoms(self) -> GeometrySocket:
        raise AttributeError(
            "Density trees have no input geometry; start each branch with "
            "`tree.import_volume()`"
        )

    def import_volume(
        self,
        name: str = "density",
        index: int = 0,
        follow_frame: bool | None = None,
        frame_start: int | None = None,
    ) -> ImportVolume:
        """
        An `Import Volume` node that reads the density's ``.vdb`` files.

        Start a branch of the tree with one, inside the tree:

        >>> with density.tree as tree:
        ...     tree.import_volume() >> DensityStyleWire() >> tree.join

        Parameters
        ----------
        name : str, optional
            Name the grid is kept under for `Color Sample Grid`, by default "density".
            Give each branch that imports a different grid its own name.
        index : int, optional
            Which of several files to import, when not following the frame.
        follow_frame : bool | None, optional
            Whether the scene frame picks the file, to play the files as an
            animation. Defaults to how the density was loaded.
        frame_start : int | None, optional
            Scene frame that shows the first file when following the frame. Defaults
            to how the density was loaded, or the scene's start frame.

        Returns
        -------
        ImportVolume
            The node, added to the tree being built.
        """
        from .grids import Grids

        density = self._entity
        assert isinstance(density, Grids), "importing a volume needs its density"
        if follow_frame is None:
            follow_frame = density._follow_frame
        if frame_start is None:
            frame_start = density._frame_start
        if frame_start is None:
            scene = bpy.context.scene
            assert scene is not None
            frame_start = scene.frame_start
        node = ImportVolume(
            path=density.vdb_path,
            index=index,
            follow_frame=follow_frame,
            frame_start=frame_start,
            count=len(density.files_vdb),
            name=name,
        )
        # saving makes the Path relative (`//`), which Import VDB inside a linked
        # group doesn't resolve (Blender 5.2, 5.3), so the group is kept local
        group = node.node.node_tree
        if group is not None and group.library is not None:
            node.node.node_tree = group.make_local()
        return node

    @contextmanager
    def reset(self, output: str = "Geometry") -> Generator[ResetSockets]:
        """
        Clear the tree and build within it.

        Yields ``(volume, join)``: the volume is the output of an `Import Volume`
        node reading the density's files as they were loaded, removed again if
        nothing is built from it. More branches can start with
        `import_volume()`.
        """
        with self:
            self.clear()
            self.outputs.geometry(output)
            volume = self.import_volume()
            yield ResetSockets(volume.o.volume, self.join)
            if not volume.node.outputs[0].links:
                self.tree.nodes.remove(volume.node)


class Density(MolecularEntity, metaclass=ABCMeta):
    """
    Abstract base class for molecular density objects.

    """

    _tree_class = DensityTree
    _entity_type = EntityType.DENSITY

    @property
    def tree(self) -> DensityTree:
        tree = super().tree
        assert isinstance(tree, DensityTree)
        return tree

    @tree.setter
    def tree(self, value) -> None:
        MolecularEntity.tree.fset(self, value)

    def __init__(self, file_path: Union[str, Path]):
        super().__init__()
        self.file_path: Path = bl.path_resolve(file_path)
        self.file_vdb: str
        self.annotations = DensityAnnotationManager(self)

    @overload
    def named_attribute(
        self, name: str, evaluate: bool = True, *, atype: None = None
    ) -> np.ndarray: ...
    @overload
    def named_attribute[T: databpy.AttributeArray](
        self, name: str, evaluate: bool = True, *, atype: type[T]
    ) -> T: ...
    def named_attribute(
        self,
        name: str,
        evaluate: bool = True,
        *,
        atype: type[databpy.AttributeArray] | None = None,
    ) -> np.ndarray:
        # volume objects can't be evaluated to a mesh data-block directly, but the
        # evaluated geometry holds the mesh their node tree generates
        return databpy.GeometrySet(self.object).named_attribute(name, atype=atype)

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
