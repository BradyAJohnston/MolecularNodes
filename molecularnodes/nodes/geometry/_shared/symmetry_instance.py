# Node group "Symmetry Instance" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    RotationSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputGeometry, InputRotation, InputVector


class SymmetryInstance(CustomGeometryGroup):
    """
    Symmetry Instance

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to place a copy of under each symmetry operator
    points : InputGeometry
        One point per symmetry operator
    rotation : InputRotation
        Rotation of each operator, evaluated on the points
    translation : InputVector
        Translation of each operator, evaluated on the points
    centre : InputVector
        Point the rotations are applied about
    factor : InputFloat
        0 places every copy back on the original, 1 builds the full symmetry
    stagger : InputFloat
        Delay each copy by its order, so at 1 the last copy only starts moving as the first finishes

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to place a copy of under each symmetry operator
    i.points : GeometrySocket
        One point per symmetry operator
    i.rotation : RotationSocket
        Rotation of each operator, evaluated on the points
    i.translation : VectorSocket
        Translation of each operator, evaluated on the points
    i.centre : VectorSocket
        Point the rotations are applied about
    i.factor : FloatSocket
        0 places every copy back on the original, 1 builds the full symmetry
    i.stagger : FloatSocket
        Delay each copy by its order, so at 1 the last copy only starts moving as the first finishes

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Symmetry Instance"
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to place a copy of under each symmetry operator"""
        points: GeometrySocket
        """One point per symmetry operator"""
        rotation: RotationSocket
        """Rotation of each operator, evaluated on the points"""
        translation: VectorSocket
        """Translation of each operator, evaluated on the points"""
        centre: VectorSocket
        """Point the rotations are applied about"""
        factor: FloatSocket
        """0 places every copy back on the original, 1 builds the full symmetry"""
        stagger: FloatSocket
        """Delay each copy by its order, so at 1 the last copy only starts moving as the first finishes"""

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
        geometry: InputGeometry = None,
        points: InputGeometry = None,
        rotation: InputRotation = None,
        translation: InputVector = None,
        centre: InputVector = None,
        factor: InputFloat = 1.0,
        stagger: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Points": points,
                "Rotation": rotation,
                "Translation": translation,
                "Centre": centre,
                "Factor": factor,
                "Stagger": stagger,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry",
            description="Geometry to place a copy of under each symmetry operator",
        )
        points = tree.inputs.geometry(
            "Points", description="One point per symmetry operator"
        )
        rotation = tree.inputs.rotation(
            "Rotation",
            (0.0, 0.0, 0.0),
            description="Rotation of each operator, evaluated on the points",
            hide_value=True,
        )
        translation = tree.inputs.vector(
            "Translation",
            (0.0, 0.0, 0.0),
            description="Translation of each operator, evaluated on the points",
            hide_value=True,
        )
        centre = tree.inputs.vector(
            "Centre",
            (0.0, 0.0, 0.0),
            description="Point the rotations are applied about",
            subtype="XYZ",
        )
        factor = tree.inputs.float(
            "Factor",
            1.0,
            description="0 places every copy back on the original, 1 builds the full symmetry",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        stagger = tree.inputs.float(
            "Stagger",
            0.0,
            description="Delay each copy by its order, so at 1 the last copy only starts moving as the first finishes",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        instances = tree.outputs.geometry("Instances")

        with g.Frame("Per-copy factor"):
            math_1 = g.Math.maximum(
                g.DomainSize(geometry=points, component="POINTCLOUD").o.point_count - 1,
                1.0,
            )
            math_2 = g.Index().o.index / math_1 * stagger
            clamp = ((factor - math_2) / (1.0 - math_2).max(0.0001)).clamp()
            _string = g.String(
                string="Each copy waits for its share of the timeline before it starts moving: with Stagger at 0 every copy shares the same Factor, at 1 the last copy only begins as the first one finishes."
            )
        with g.Frame("Partial operator"):
            rotation_to_axis_angle = rotation.to_axis_angle()
            axis_angle_to_rotation = g.AxisAngleToRotation(
                axis=rotation_to_axis_angle.axis,
                angle=rotation_to_axis_angle.angle * clamp,
            )
            _string_1 = g.String(
                string="Rotating about the operator's own axis by a fraction of its angle is the exact interpolation from the identity to that operator. A rotation about a centre c is R p + (c - R c), so the centre is folded into the position; at factor 0 that is zero and every copy sits on the original."
            )
            vector_math = (
                centre - centre.rotate(axis_angle_to_rotation) + translation * clamp
            )
        (
            points
            >> g.SetPosition(position=vector_math)
            >> g.InstanceOnPoints(instance=geometry, rotation=axis_angle_to_rotation)
            >> g.StoreNamedAttribute.instance.integer(name="sym_id", value=g.Index())
            >> instances
        )
