# Node-group asset 'Curve Visualize' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
import bpy
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry, InputVector
from ._shared.mn_units import MNUnits
from ._shared.unit_convert import UnitConvert
from .curve_rotation import CurveRotation
from .primitive_arrow import PrimitiveArrow
from .primitive_gimbal import PrimitiveGimbal
from .set_color import SetColor


class CurveVisualize(AssetGeometryGroup):
    """
    Curve Visualize

    Parameters
    ----------
    curve : InputGeometry
        Curve
    selection : InputBoolean
        Selection
    position : InputVector
        Position
    normal : InputVector
        Normal
    handles : InputBoolean
        Handles
    value : InputFloat
        Value
    arrow_size : InputFloat
        Arrow Size

    Inputs
    ------
    i.curve : GeometrySocket
        Curve
    i.selection : BooleanSocket
        Selection
    i.position : VectorSocket
        Position
    i.normal : VectorSocket
        Normal
    i.handles : BooleanSocket
        Handles
    i.value : FloatSocket
        Value
    i.arrow_size : FloatSocket
        Arrow Size

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Curve Visualize"
    _asset_name = "Curve Visualize"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.curve_visualize"}

    class _Inputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""
        selection: BooleanSocket
        """Selection"""
        position: VectorSocket
        """Position"""
        normal: VectorSocket
        """Normal"""
        handles: BooleanSocket
        """Handles"""
        value: FloatSocket
        """Value"""
        arrow_size: FloatSocket
        """Arrow Size"""

    class _Outputs(SocketAccessor):
        instances: GeometrySocket
        """Instances"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        curve: InputGeometry = None,
        selection: InputBoolean = True,
        position: InputVector = None,
        normal: InputVector = None,
        handles: InputBoolean = False,
        value: InputFloat = 3.0,
        arrow_size: InputFloat = 2.0,
    ):
        super().__init__(
            **{
                "Curve": curve,
                "Selection": selection,
                "Position": position,
                "Normal": normal,
                "Handles": handles,
                "Value": value,
                "Arrow Size": arrow_size,
            }
        )

    def _build_group(self, tree):
        curve = tree.inputs.geometry("Curve")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        position = tree.inputs.vector(
            "Position", (0.0, 0.0, 0.0), hide_value=True, default_input="POSITION"
        )
        normal = tree.inputs.vector(
            "Normal", (0.0, 0.0, 1.0), hide_value=True, default_input="NORMAL"
        )
        handles = tree.inputs.boolean("Handles", False)
        _value = tree.inputs.float(
            "Value", 3.0, min_value=-10_000.0, max_value=10_000.0
        )
        arrow_size = tree.inputs.float(
            "Arrow Size", 2.0, min_value=-10_000.0, max_value=10_000.0
        )
        instances = tree.outputs.geometry("Instances")

        curve_handle_positions = g.CurveHandlePositions(relative=True)
        group = PrimitiveArrow(
            vertices=3,
            value=(0.1056426, 0.800023, 0.7937081, 1.0),
            material=bpy.data.materials["MN Ambient Occlusion"],
        )
        capture = g.CaptureAttribute.point(geometry=curve)
        selection_1 = capture.items.boolean("Selection", selection)
        position_1 = capture.items.vector("Position", position)
        normal_1 = capture.items.vector("Normal", normal)
        group_1 = UnitConvert(from_=2.0)
        _group_2 = CurveRotation(normal=normal_1.output)
        group_3 = PrimitiveArrow(
            vertices=3,
            value=(0.8000315, 0.49981865, 0.01965215, 1.0),
            material=bpy.data.materials["MN Ambient Occlusion"],
        )
        group_4 = UnitConvert(from_=2.0)
        set_spline_type = g.SetSplineType.bezier(
            g.SeparateComponents(geometry=capture.o.geometry).o.curve
        )
        set_spline_type.node.mute = True
        set_handle_type = g.SetHandleType(curve=set_spline_type)
        set_handle_type.node.mute = True
        group_5 = SetColor(atoms=set_handle_type, color=(0.0, 0.0, 0.0, 1.0))
        set_position = g.SetPosition(
            geometry=group_5, selection=selection_1.output, position=position_1.output
        )
        instance_on_points = g.InstanceOnPoints(
            points=set_position,
            instance=PrimitiveGimbal(
                vertices=3, material=bpy.data.materials["MN Ambient Occlusion"]
            ),
            rotation=g.NamedAttribute.quaternion("rotation").o.attribute,
            scale=MNUnits(value=arrow_size).o.angstrom,
        )
        instance_on_points_1 = g.InstanceOnPoints(
            points=group_5,
            instance=group,
            rotation=g.AlignRotationToVector(vector=curve_handle_positions.o.left),
            scale=g.CombineXYZ(
                x=group_1, y=group_1, z=curve_handle_positions.o.left.length()
            ),
        )
        instance_on_points_2 = g.InstanceOnPoints(
            points=group_5,
            instance=group_3,
            rotation=g.AlignRotationToVector(vector=curve_handle_positions.o.right),
            scale=g.CombineXYZ(
                x=group_4, y=group_4, z=curve_handle_positions.o.right.length()
            ),
        )
        join_geometry = g.JoinGeometry(
            geometry=(
                g.JoinGeometry(geometry=(instance_on_points_1, instance_on_points_2)),
                instance_on_points_1,
            )
        )
        join_geometry_1 = g.JoinGeometry(
            geometry=(
                handles.switch.geometry(true=join_geometry),
                instance_on_points,
                capture.o.geometry,
            )
        )

        join_geometry_1 >> instances


ASSET = CurveVisualize

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}

DATABLOCK_DEPENDENCIES = {
    "materials": ("MN Ambient Occlusion",),
}
