# Node-group asset "Symmetry Dihedral" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputGeometry, InputInteger, InputVector
from ._shared.symmetry_instance import SymmetryInstance


class SymmetryDihedral(AssetGeometryGroup):
    """
    Build a dihedral (Dn) assembly: a ring of n copies plus the same ring flipped by a perpendicular two-fold, giving 2n copies

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to replicate around the axis
    order : InputInteger
        Number of copies in each ring (the n of Dn), giving 2n copies in total
    axis : InputVector
        Direction of the main symmetry axis
    centre : InputVector
        Point the symmetry axes pass through
    offset : InputFloat
        Rotation of the flipped ring about the axis, which sets where the perpendicular two-fold lies
    factor : InputFloat
        0 places every copy back on the original, 1 builds the full symmetry
    stagger : InputFloat
        Delay each copy by its order, so at 1 the last copy only starts moving as the first finishes

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to replicate around the axis
    i.order : IntegerSocket
        Number of copies in each ring (the n of Dn), giving 2n copies in total
    i.axis : VectorSocket
        Direction of the main symmetry axis
    i.centre : VectorSocket
        Point the symmetry axes pass through
    i.offset : FloatSocket
        Rotation of the flipped ring about the axis, which sets where the perpendicular two-fold lies
    i.factor : FloatSocket
        0 places every copy back on the original, 1 builds the full symmetry
    i.stagger : FloatSocket
        Delay each copy by its order, so at 1 the last copy only starts moving as the first finishes

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Symmetry Dihedral"
    _asset_name = "Symmetry Dihedral"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Build a dihedral (Dn) assembly: a ring of n copies plus the same ring flipped by a perpendicular two-fold, giving 2n copies"
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to replicate around the axis"""
        order: IntegerSocket
        """Number of copies in each ring (the n of Dn), giving 2n copies in total"""
        axis: VectorSocket
        """Direction of the main symmetry axis"""
        centre: VectorSocket
        """Point the symmetry axes pass through"""
        offset: FloatSocket
        """Rotation of the flipped ring about the axis, which sets where the perpendicular two-fold lies"""
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
        order: InputInteger = 3,
        axis: InputVector = None,
        centre: InputVector = None,
        offset: InputFloat = 0.0,
        factor: InputFloat = 1.0,
        stagger: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Order": order,
                "Axis": axis,
                "Centre": centre,
                "Offset": offset,
                "Factor": factor,
                "Stagger": stagger,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to replicate around the axis"
        )
        order = tree.inputs.integer(
            "Order",
            3,
            description="Number of copies in each ring (the n of Dn), giving 2n copies in total",
            min_value=1,
            max_value=1000,
        )
        axis = tree.inputs.vector(
            "Axis",
            (0.0, 0.0, 1.0),
            description="Direction of the main symmetry axis",
            subtype="XYZ",
        )
        centre = tree.inputs.vector(
            "Centre",
            (0.0, 0.0, 0.0),
            description="Point the symmetry axes pass through",
            subtype="XYZ",
        )
        offset = tree.inputs.float(
            "Offset",
            0.0,
            description="Rotation of the flipped ring about the axis, which sets where the perpendicular two-fold lies",
            subtype="ANGLE",
        )
        with tree.inputs.panel("Animate", default_closed=True):
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

        index = g.Index()
        with g.Frame("Dn operators"):
            axis_angle_to_rotation = g.AxisAngleToRotation(
                axis=axis,
                angle=index.o.index.modulo(order) * (math.tau / order) + offset,
            )
            axis_angle_to_rotation_1 = g.AxisAngleToRotation(
                axis=g.RotateVector(
                    rotation=g.AlignRotationToVector(vector=axis),
                    vector=(0.0, 1.0, 0.0),
                ),
                angle=math.pi,
            )
            rotate_rotation = axis_angle_to_rotation.o.rotation.rotate(
                axis_angle_to_rotation_1, rotation_space="LOCAL"
            )
            _string = g.String(
                string="Points 0..n-1 are the first ring, n..2n-1 the same ring flipped by a two-fold perpendicular to the axis. The two-fold is Y carried onto the plane perpendicular to the axis; Offset rotates the flipped ring about the axis, which moves that two-fold by half the offset."
            )
            switch = (index >= order).switch.rotation(
                g.AxisAngleToRotation(
                    axis=axis, angle=index.o.index.modulo(order) * (math.tau / order)
                ),
                rotate_rotation,
            )
        (
            SymmetryInstance(
                geometry=geometry,
                points=g.Points(count=order * 2, radius=0.1),
                rotation=switch,
                centre=centre,
                factor=factor,
                stagger=stagger,
            )
            >> instances
        )


ASSET = SymmetryDihedral

ASSET_METADATA = {
    "description": "Build a dihedral (Dn) assembly: a ring of n copies plus the same ring flipped by a perpendicular two-fold, giving 2n copies",
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
