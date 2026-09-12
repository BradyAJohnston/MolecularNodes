# Node-group asset "Select Sphere" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .boolean_andor import BooleanAndOr


class SelectSphere(AssetGeometryGroup):
    """
    Select Sphere

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection
    object : InputObject
        The position from this Object will be used for the distance calculation for the selection. By default an Empty object is used, but any object can be used in principle

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection
    i.object : ObjectSocket
        The position from this Object will be used for the distance calculation for the selection. By default an Empty object is used, but any object can be used in principle

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = "Select Sphere"
    _asset_name = "Select Sphere"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_sphere"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""
        object: ObjectSocket
        """The position from this Object will be used for the distance calculation for the selection. By default an Empty object is used, but any object can be used in principle"""

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
            bpy.data.objects.get("select_sphere"),
            description="The position from this Object will be used for the distance calculation for the selection. By default an Empty object is used, but any object can be used in principle",
            optional_label=True,
        )
        selection = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )
        inverted = tree.outputs.boolean(
            "Inverted", description="The inverse of the calculated selection"
        )

        object_info = g.ObjectInfo(
            object=object, as_instance=True, transform_space="RELATIVE"
        )
        compare = g.Compare.float.less_than(
            g.Position().o.position.distance(object_info.o.location),
            abs(object_info.o.scale),
        )
        group = BooleanAndOr(and_=and_, or_=or_, boolean=compare)

        group >> selection
        group.o.inverted >> inverted


ASSET = SelectSphere

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}

DATABLOCK_DEPENDENCIES = {
    "objects": ("select_sphere",),
}
