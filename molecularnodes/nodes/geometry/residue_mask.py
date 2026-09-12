# Node-group asset "Residue Mask" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputInteger
from .atom_name import AtomName
from .group_pick_vector import GroupPickVector
from .ures_id import UResID


class ResidueMask(AssetGeometryGroup):
    """
    Residue Mask

    Parameters
    ----------
    atom_name : InputInteger
        Atom to pick from the group

    Inputs
    ------
    i.atom_name : IntegerSocket
        Atom to pick from the group

    Outputs
    -------
    o.is_valid : BooleanSocket
        Group contains only one occurrance of the selected atom. None or more than one returns False
    o.index : IntegerSocket
        Index for the group's atom with specified name, returns -1 if not valid
    o.position : VectorSocket
        Position of the picked point in the group, returns (0, 0, 0) if not valid
    """

    _name = "Residue Mask"
    _asset_name = "Residue Mask"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.residue_mask"}

    class _Inputs(SocketAccessor):
        atom_name: IntegerSocket
        """Atom to pick from the group"""

    class _Outputs(SocketAccessor):
        is_valid: BooleanSocket
        """Group contains only one occurrance of the selected atom. None or more than one returns False"""
        index: IntegerSocket
        """Index for the group's atom with specified name, returns -1 if not valid"""
        position: VectorSocket
        """Position of the picked point in the group, returns (0, 0, 0) if not valid"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atom_name: InputInteger = 1,
    ):
        super().__init__(**{"atom_name": atom_name})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atom_name = tree.inputs.integer(
            "atom_name", 1, description="Atom to pick from the group", min_value=2
        )
        is_valid = tree.outputs.boolean(
            "Is Valid",
            description="Group contains only one occurrance of the selected atom. None or more than one returns False",
        )
        index = tree.outputs.integer(
            "Index",
            description="Index for the group's atom with specified name, returns -1 if not valid",
        )
        position = tree.outputs.vector(
            "Position",
            description="Position of the picked point in the group, returns (0, 0, 0) if not valid",
        )

        group = GroupPickVector(
            pick=g.Compare.integer.equal(AtomName(), atom_name),
            group_id=UResID().o.ures_id,
        )

        group >> is_valid
        group.o.index >> index
        group.o.vector >> position


ASSET = ResidueMask

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
