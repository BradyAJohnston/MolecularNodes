# Node-group asset "Select Cube" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
import bpy
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ObjectSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputObject
from .between_vector import BetweenVector
from .boolean_andor import BooleanAndOr


class SelectCube(AssetGeometryGroup):
    """
    Select Cube

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection
    object : InputObject
        The position, rotation and scale from this Object will be used for the distance calculation for the selection. By default an Empty object is used, but any object can be used in principle

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection
    i.object : ObjectSocket
        The position, rotation and scale from this Object will be used for the distance calculation for the selection. By default an Empty object is used, but any object can be used in principle

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = "Select Cube"
    _asset_name = "Select Cube"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_cube"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""
        object: ObjectSocket
        """The position, rotation and scale from this Object will be used for the distance calculation for the selection. By default an Empty object is used, but any object can be used in principle"""

    class _Outputs(SocketAccessor):
        selection: BooleanSocket
        """The calculated selection"""
        inverted: BooleanSocket
        """The inverse of the calculated selection"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        and_: InputBoolean = True,
        or_: InputBoolean = False,
        object: InputObject = None,
    ):
        super().__init__(**{"And": and_, "Or": or_, "Object": object})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        and_ = tree.inputs.boolean(
            "And",
            True,
            description="The resulting selection must overlap with this input selection",
            hide_value=True,
        )
        or_ = tree.inputs.boolean(
            "Or",
            False,
            description="The resulting selection can be calculated from this node or be from this input selection",
            hide_value=True,
        )
        object = tree.inputs.object(
            "Object",
            bpy.data.objects.get("select_cube"),
            description="The position, rotation and scale from this Object will be used for the distance calculation for the selection. By default an Empty object is used, but any object can be used in principle",
            optional_label=True,
        )
        selection = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )
        inverted = tree.outputs.boolean(
            "Inverted", description="The inverse of the calculated selection"
        )

        invert_matrix = g.ObjectInfo(
            object=object, as_instance=True, transform_space="RELATIVE"
        ).o.transform.invert()
        group = BetweenVector(
            value=g.ProjectPoint(vector=g.Position(), transform=invert_matrix),
            lower=(-1.0, -1.0, -1.0),
            upper=(1.0, 1.0, 1.0),
        )
        group_1 = BooleanAndOr(and_=and_, or_=or_, boolean=group)

        group_1 >> selection
        group_1.o.inverted >> inverted


ASSET = SelectCube

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}

DATABLOCK_DEPENDENCIES = {
    "objects": ("select_cube",),
}
