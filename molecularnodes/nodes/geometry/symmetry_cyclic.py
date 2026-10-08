# Node-group asset "Symmetry Cyclic" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
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


class SymmetryCyclic(AssetGeometryGroup):
    """
    Build a cyclic (Cn) assembly: n copies evenly spaced around one axis

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to replicate around the axis
    order : InputInteger
        Number of copies evenly spaced around the axis (the n of Cn)
    axis : InputVector
        Direction of the symmetry axis
    centre : InputVector
        Point the symmetry axis passes through
    animate : InputFloat
        0 places a copy back on the original, 1 builds the full symmetry. Evaluated on each copy, so a field such as Stagger Value moves the copies one after another

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to replicate around the axis
    i.order : IntegerSocket
        Number of copies evenly spaced around the axis (the n of Cn)
    i.axis : VectorSocket
        Direction of the symmetry axis
    i.centre : VectorSocket
        Point the symmetry axis passes through
    i.animate : FloatSocket
        0 places a copy back on the original, 1 builds the full symmetry. Evaluated on each copy, so a field such as Stagger Value moves the copies one after another

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Symmetry Cyclic"
    _asset_name = "Symmetry Cyclic"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Build a cyclic (Cn) assembly: n copies evenly spaced around one axis"
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to replicate around the axis"""
        order: IntegerSocket
        """Number of copies evenly spaced around the axis (the n of Cn)"""
        axis: VectorSocket
        """Direction of the symmetry axis"""
        centre: VectorSocket
        """Point the symmetry axis passes through"""
        animate: FloatSocket
        """0 places a copy back on the original, 1 builds the full symmetry. Evaluated on each copy, so a field such as Stagger Value moves the copies one after another"""

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
        animate: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Order": order,
                "Axis": axis,
                "Centre": centre,
                "Animate": animate,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to replicate around the axis"
        )
        order = tree.inputs.integer(
            "Order",
            3,
            description="Number of copies evenly spaced around the axis (the n of Cn)",
            min_value=1,
            max_value=1000,
        )
        axis = tree.inputs.vector(
            "Axis",
            (0.0, 0.0, 1.0),
            description="Direction of the symmetry axis",
            subtype="XYZ",
        )
        centre = tree.inputs.vector(
            "Centre",
            (0.0, 0.0, 0.0),
            description="Point the symmetry axis passes through",
            subtype="XYZ",
        )
        animate = tree.inputs.float(
            "Animate",
            1.0,
            description="0 places a copy back on the original, 1 builds the full symmetry. Evaluated on each copy, so a field such as Stagger Value moves the copies one after another",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        instances = tree.outputs.geometry("Instances")

        with g.Frame("Cn operators"):
            axis_angle_to_rotation = g.AxisAngleToRotation(
                axis=axis, angle=g.Index().o.index * (math.tau / order)
            )
            _string = g.String(
                string="Copy k is rotated k * 360 / n degrees about the axis. The rotation is passed to Symmetry Instance as a field and evaluated on the points."
            )
        (
            SymmetryInstance(
                geometry=geometry,
                points=g.Points(count=order),
                rotation=axis_angle_to_rotation,
                centre=centre,
                animate=animate,
            )
            >> instances
        )


ASSET = SymmetryCyclic

ASSET_METADATA = {
    "description": "Build a cyclic (Cn) assembly: n copies evenly spaced around one axis",
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
