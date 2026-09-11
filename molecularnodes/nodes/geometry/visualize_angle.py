# Node-group asset "Visualize Angle" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
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


class VisualizeAngle(AssetGeometryGroup):
    """
    Visualize Angle

    Parameters
    ----------
    points : InputGeometry
        Points
    selection : InputBoolean
        Selection
    position : InputVector
        Position
    angle : InputFloat
        Angle
    length : InputFloat
        Length
    up : InputVector
        Up
    axis : InputVector
        Axis
    radius : InputFloat
        Radius

    Inputs
    ------
    i.points : GeometrySocket
        Points
    i.selection : BooleanSocket
        Selection
    i.position : VectorSocket
        Position
    i.angle : FloatSocket
        Angle
    i.length : FloatSocket
        Length
    i.up : VectorSocket
        Up
    i.axis : VectorSocket
        Axis
    i.radius : FloatSocket
        Radius

    Outputs
    -------
    o.mesh : GeometrySocket
        Mesh
    o.curve : GeometrySocket
        Curve
    """

    _name = "Visualize Angle"
    _asset_name = "Visualize Angle"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        points: GeometrySocket
        """Points"""
        selection: BooleanSocket
        """Selection"""
        position: VectorSocket
        """Position"""
        angle: FloatSocket
        """Angle"""
        length: FloatSocket
        """Length"""
        up: VectorSocket
        """Up"""
        axis: VectorSocket
        """Axis"""
        radius: FloatSocket
        """Radius"""

    class _Outputs(SocketAccessor):
        mesh: GeometrySocket
        """Mesh"""
        curve: GeometrySocket
        """Curve"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        points: InputGeometry = None,
        selection: InputBoolean = True,
        position: InputVector = None,
        angle: InputFloat = 0.5,
        length: InputFloat = 0.005,
        up: InputVector = None,
        axis: InputVector = None,
        radius: InputFloat = 0.2,
    ):
        super().__init__(
            **{
                "Points": points,
                "Selection": selection,
                "Position": position,
                "Angle": angle,
                "Length": length,
                "Up": up,
                "Axis": axis,
                "Radius": radius,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        points = tree.inputs.geometry("Points")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        position = tree.inputs.vector(
            "Position", (0.0, 0.0, 0.0), hide_value=True, default_input="POSITION"
        )
        angle = tree.inputs.float("Angle", 0.5, min_value=-10_000.0, max_value=10_000.0)
        length = tree.inputs.float(
            "Length", 0.005, min_value=-10_000.0, max_value=10_000.0
        )
        up = tree.inputs.vector(
            "Up", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0
        )
        axis = tree.inputs.vector("Axis", (0.0, 0.0, 1.0))
        radius = tree.inputs.float("Radius", 0.2, min_value=0.0, subtype="DISTANCE")
        mesh = tree.outputs.geometry("Mesh")
        curve = tree.outputs.geometry("Curve")

        named_attribute = g.NamedAttribute.float("radius")
        resample_curve = g.ResampleCurve(
            curve=g.CurveLine(end=(0.0, 0.0, 0.0)), length=0.1, keep_last_segment=True
        )
        store_named_attribute = (
            points
            >> g.InstanceOnPoints(selection=selection, instance=resample_curve)
            >> g.RealizeInstances(realize_to_point_domain=True)
            >> g.StoreNamedAttribute.point.float(
                name="angle", value=angle * g.SplineParameter().o.factor * -1.0
            )
        )
        rotate_vector = (up.normalize() * length).rotate(
            g.AxisAngleToRotation(
                axis=axis, angle=g.NamedAttribute.float("angle").o.attribute
            )
        )
        capture = g.CaptureAttribute.point(geometry=store_named_attribute)
        vector = capture.items.vector("Vector", rotate_vector)
        set_curve_radius = (
            capture.o.geometry
            >> g.SetCurveNormal(normal=vector.output, mode="Free")
            >> g.SetPosition(position=position, offset=vector.output)
            >> g.SetCurveRadius(radius=radius)
        )
        curve_to_mesh = g.CurveToMesh(
            curve=set_curve_radius,
            profile_curve=g.CurveCircle(resolution=6, radius=0.01),
            scale=named_attribute.o.exists.switch.float(
                1.0, named_attribute.o.attribute
            ),
            fill_caps=True,
        )

        curve_to_mesh >> mesh
        set_curve_radius >> curve


ASSET = VisualizeAngle

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
