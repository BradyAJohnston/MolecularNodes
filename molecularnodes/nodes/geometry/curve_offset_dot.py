# Node-group asset 'Curve Offset Dot' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputInteger, InputMenu, InputVector
from .offset_vector import OffsetVector


class CurveOffsetDot(AssetGeometryGroup):
    """
    Curve Offset Dot

    Parameters
    ----------
    normal : InputVector
        Normal
    offset : InputInteger
        Offset
    threshold_direction : InputMenu | Literal["Less Than", "Greater Than"]
        Threshold Direction
    threshold_cutoff : InputFloat
        Threshold Cutoff
    rotation_axis : InputVector
        Rotation Axis
    rotation_amount : InputFloat
        Rotation Amount

    Inputs
    ------
    i.normal : VectorSocket
        Normal
    i.offset : IntegerSocket
        Offset
    i.threshold_direction : MenuSocket
        Threshold Direction
    i.threshold_cutoff : FloatSocket
        Threshold Cutoff
    i.rotation_axis : VectorSocket
        Rotation Axis
    i.rotation_amount : FloatSocket
        Rotation Amount

    Outputs
    -------
    o.thresholded : BooleanSocket
        Thresholded
    o.leading : IntegerSocket
        Leading
    o.rotation : RotationSocket
        Rotation
    """

    _name = "Curve Offset Dot"
    _asset_name = "Curve Offset Dot"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.curve_offset_dot"}

    class _Inputs(SocketAccessor):
        normal: VectorSocket
        """Normal"""
        offset: IntegerSocket
        """Offset"""
        threshold_direction: MenuSocket
        """Threshold Direction"""
        threshold_cutoff: FloatSocket
        """Threshold Cutoff"""
        rotation_axis: VectorSocket
        """Rotation Axis"""
        rotation_amount: FloatSocket
        """Rotation Amount"""

    class _Outputs(SocketAccessor):
        thresholded: BooleanSocket
        """Thresholded"""
        leading: IntegerSocket
        """Leading"""
        rotation: RotationSocket
        """Rotation"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        normal: InputVector = None,
        offset: InputInteger = -1,
        threshold_direction: InputMenu
        | Literal["Less Than", "Greater Than"] = "Less Than",
        threshold_cutoff: InputFloat = -0.9,
        rotation_axis: InputVector = None,
        rotation_amount: InputFloat = math.pi,
    ):
        super().__init__(
            **{
                "Normal": normal,
                "Offset": offset,
                "Threshold Direction": threshold_direction,
                "Threshold Cutoff": threshold_cutoff,
                "Rotation Axis": rotation_axis,
                "Rotation Amount": rotation_amount,
            }
        )

    def _build_group(self, tree):
        normal = tree.inputs.vector("Normal", (0.0, 0.0, 0.0), default_input="NORMAL")
        offset = tree.inputs.integer("Offset", -1, min_value=-2147483647)
        with tree.inputs.panel("Threshold"):
            threshold_direction = tree.inputs.menu(
                "Threshold Direction", optional_label=True
            )
            threshold_cutoff = tree.inputs.float(
                "Threshold Cutoff", -0.9, min_value=-10_000.0, max_value=10_000.0
            )
        with tree.inputs.panel("Rotation", default_closed=True):
            rotation_axis = tree.inputs.vector("Rotation Axis", (0.0, 0.0, 1.0))
            rotation_amount = tree.inputs.float(
                "Rotation Amount", math.pi, min_value=-10_000.0, max_value=10_000.0
            )
        thresholded = tree.outputs.boolean("Thresholded")
        leading = tree.outputs.integer("Leading")
        rotation = tree.outputs.rotation("Rotation")

        vector_math = normal.dot(OffsetVector(vector=normal, offset=offset))
        menu_switch = g.MenuSwitch.boolean(
            threshold_direction,
            {
                "Less Than": vector_math < threshold_cutoff,
                "Greater Than": vector_math > threshold_cutoff,
            },
        )
        accumulate_field = g.AccumulateField.point.integer(
            menu_switch.o.output, g.CurveOfPoint().o.curve_index
        )
        axis_angle_to_rotation = g.AxisAngleToRotation(
            axis=rotation_axis, angle=accumulate_field.o.trailing * rotation_amount
        )

        menu_switch >> thresholded
        accumulate_field.o.trailing >> leading
        axis_angle_to_rotation >> rotation

        threshold_direction.default_value = "Less Than"


ASSET = CurveOffsetDot

ASSET_METADATA = {
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
