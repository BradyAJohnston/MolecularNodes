from abc import ABCMeta
from contextlib import contextmanager
from enum import Enum
from typing import Iterator, List, NamedTuple, cast
import bpy
from bpy.types import GeometryNodeTree
from databpy import (
    BlenderObject,
)
from nodebpy import geometry as g
from nodebpy.builder import GeometrySocket, TreeBuilder
from nodebpy.nodes.geometry import JoinGeometry
from ..blender import utils as blender_utils
from .utilities import BoolObjectMNProperty


# create a EntityType enum
# These should match the entity_type EnumProperty in props.py
class EntityType(Enum):
    MD = "md"
    MD_OXDNA = "md-oxdna"
    MD_STREAMING = "md-streaming"
    MOLECULE = "molecule"
    ENSEMBLE = "ensemble"
    DENSITY = "density"
    ENSEMBLE_STAR = "ensemble-star"
    ENSEMBLE_CELLPACK = "ensemble-cellpack"


class ResetSockets(NamedTuple):
    """The sockets to build between, as returned by `MolecularTree.reset()`."""

    atoms: GeometrySocket
    join: GeometrySocket


class MolecularTree(TreeBuilder[GeometryNodeTree]):
    def __init__(self, entity: "MolecularEntity", tree: TreeBuilder | GeometryNodeTree | str) -> None:
        self._entity = entity
        if isinstance(tree, TreeBuilder):
            super().__init__(cast(GeometryNodeTree, tree.tree))
        elif isinstance(tree, GeometryNodeTree):
            super().__init__(tree)
        else:
            super().__init__(tree)


    def _wrap(self, socket: bpy.types.NodeSocket) -> GeometrySocket:
        """Wrap an existing Blender socket as a socket bound to this tree."""
        wrapped = GeometrySocket(socket)
        wrapped._tree = self
        return wrapped

    def _interface_geometry(
        self, in_out: str
    ) -> bpy.types.NodeTreeInterfaceSocket | None:
        """The tree's first geometry interface socket in the given direction, if any."""
        for item in self.tree.interface.items_tree:
            if (
                item.item_type == "SOCKET"
                and item.in_out == in_out
                and item.socket_type == "NodeSocketGeometry"
            ):
                return item
        return None

    @property
    def atoms(self) -> GeometrySocket:
        """
        The geometry input to build branches from, adding the input if it is missing.

        >>> with mol.tree as tree:
        ...     tree.atoms >> StyleCartoon() >> tree.join
        """
        item = self._interface_geometry("INPUT")
        if item is None:
            return self.inputs.geometry("Atoms")
        return self._wrap(self._input_node().outputs[item.identifier])

    @property
    def geometry(self) -> GeometrySocket:
        """The geometry output of the tree, adding the output if it is missing."""
        item = self._interface_geometry("OUTPUT")
        if item is None:
            return self.outputs.geometry("Geometry")
        return self._wrap(self._output_node().inputs[item.identifier])

    @property
    def join(self) -> GeometrySocket:
        """
        The join that feeds the tree output, which every branch ends in.

        Adds the join if it is missing, keeping whatever already fed the output.
        """
        output = self.geometry
        links = output.socket.links
        if links:
            node = links[0].from_socket.node
            if node.bl_idname == "GeometryNodeJoinGeometry":
                return self._wrap(node.inputs[0])

        join = g.JoinGeometry()
        if links and links[0].from_socket.node.bl_idname != "NodeGroupInput":
            # real work already feeds the output, so join it rather than orphan it.
            # a bare input -> output passthrough is just the tree's default state and
            # is dropped, otherwise unstyled geometry is rendered alongside every branch
            self._wrap(links[0].from_socket) >> join
        join >> output
        return join.i.geometry

    @contextmanager
    def reset(
        self, input: str = "Atoms", output: str = "Geometry"
    ) -> Iterator[ResetSockets]:
        """
        Clear the tree back to a default state and build within it.

        Discards the existing tree, so use `with entity.tree` instead to add to it.

        >>> with mol.tree.reset() as (atoms, join):
        ...     atoms >> StyleCartoon() >> join
        """
        with self:
            self.clear()
            self.inputs.geometry(input)
            self.outputs.geometry(output)
            yield ResetSockets(self.atoms, self.join)

    def clear(self) -> None:
        self.tree.nodes.clear()
        assert self.tree.interface
        self.tree.interface.clear()


class MolecularEntity(
    BlenderObject,
    metaclass=ABCMeta,
):
    update_with_scene = BoolObjectMNProperty("update_with_scene")
    _entity_type: EntityType | None = None

    def __init__(self) -> None:
        super().__init__(obj=None)
        self._register_with_session()
        self._world_scale = 0.1
        self._tree = None

    def __getstate__(self):
        """Custom serialization."""
        state = self.__dict__.copy()
        # the cached tree wraps live Blender data which can't be pickled, and is
        # rebuilt on the next access of the `tree` property
        state["_tree"] = None
        return state

    @property
    def node_group(self) -> bpy.types.GeometryNodeTree:
        mod = self.modifier
        mod.node_group = self.tree.tree
        return mod.node_group

    @property
    def modifier(self) -> bpy.types.NodesModifier:
        for mod in self.object.modifiers:
            if mod.type == "NODES" and mod.name == "Molecular Nodes":
                return mod

        mod = self.object.modifiers.new("GeometryNodes", "NODES")
        return cast(bpy.types.NodesModifier, mod)

    @property
    def tree(self) -> MolecularTree:
        if self._tree is None:
            self._tree = MolecularTree(self, self.modifier_node_tree)
        return self._tree

    @property
    def modifier_node_tree(self) -> bpy.types.GeometryNodeTree:
        mod: bpy.types.NodesModifier = self.object.modifiers.get("Molecular Nodes")

        if mod is None:
            mod = self.object.modifiers.new("Molecular Nodes", "NODES")

        if mod.node_group is None:
            mod.node_group = g.tree().tree

        mod.node_group.is_modifier = True
        return cast(GeometryNodeTree, mod.node_group)

    def _register_with_session(self) -> None:
        bpy.context.scene.MNSession.register_entity(self)  # type: ignore

    def set_frame(self, frame: int) -> None:
        """
        Update the underlying data to correspond to changes in the scene frame.

        This method should be implemented by subclasses to update the entity's
        data based on the given frame number. This can include updating positions
        and performing other necessary calculations.

        Args:
            frame (int): The frame number to update the entity's data to.

        Raises:
            NotImplementedError: If the method is not implemented by a subclass.
        """
        raise NotImplementedError("Subclasses must implement this method")

    def _setup_modifiers(self):
        """
        Create the modifiers for the molecule.
        """
        name = f"MN_{self.name}"
        with TreeBuilder.geometry(name) as tree:
            (
                tree.inputs.geometry("Atoms")
                >> tree.outputs.geometry("Geometry")
            )
        self.object.modifiers.new("Molecular Nodes", "NODES")
        self.object.modifiers[-1].node_group = tree.tree

    def get_view(self) -> List[tuple]:
        """
        Get the 3D bounding box of the entity object

        """
        return blender_utils.get_bounding_box(self.object)

    def _get_annotation_entity_type(self) -> str:
        """
        Internal: Get entity type string for annotations.

        Subentities that derive from other entities can set this
        to the parent entity type to re-use annotations of the parent.
        By default, this returns the string value of the Entity type.

        Eg: OXDNA and StreamingTrajectory that derive from Trajectory entity
        can return EntityType.MD.value to re-use the annotations from
        the parent Trajectory entity.

        """
        return self._entity_type.value
