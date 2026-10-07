# Node-group asset "Symmetry Helical" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .angstrom_to_world import AngstromToWorld


class SymmetryHelical(AssetGeometryGroup):
    """
    Build a helical assembly: each subunit advances along the axis by the rise and rotates about it by the twist

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to replicate along the helix
    count : InputInteger
        Number of subunits to place along the helix, each repeated Axial Symmetry times around the axis
    rise : InputFloat
        Distance advanced along the axis per subunit, in Angstrom
    twist : InputFloat
        Angle rotated about the axis per subunit
    axial_symmetry : InputInteger
        Cn rotational symmetry about the helical axis: each subunit is repeated this many times evenly around the axis, as in filaments built from several protofilaments
    axis : InputVector
        Direction of the helical axis
    centre : InputVector
        Point the helical axis passes through
    animate : InputFloat
        0 places a copy back on the original, 1 builds the full symmetry. Evaluated on each copy, so a field such as Stagger Value moves the copies one after another

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to replicate along the helix
    i.count : IntegerSocket
        Number of subunits to place along the helix, each repeated Axial Symmetry times around the axis
    i.rise : FloatSocket
        Distance advanced along the axis per subunit, in Angstrom
    i.twist : FloatSocket
        Angle rotated about the axis per subunit
    i.axial_symmetry : IntegerSocket
        Cn rotational symmetry about the helical axis: each subunit is repeated this many times evenly around the axis, as in filaments built from several protofilaments
    i.axis : VectorSocket
        Direction of the helical axis
    i.centre : VectorSocket
        Point the helical axis passes through
    i.animate : FloatSocket
        0 places a copy back on the original, 1 builds the full symmetry. Evaluated on each copy, so a field such as Stagger Value moves the copies one after another

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Symmetry Helical"
    _asset_name = "Symmetry Helical"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Build a helical assembly: each subunit advances along the axis by the rise and rotates about it by the twist"
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to replicate along the helix"""
        count: IntegerSocket
        """Number of subunits to place along the helix, each repeated Axial Symmetry times around the axis"""
        rise: FloatSocket
        """Distance advanced along the axis per subunit, in Angstrom"""
        twist: FloatSocket
        """Angle rotated about the axis per subunit"""
        axial_symmetry: IntegerSocket
        """Cn rotational symmetry about the helical axis: each subunit is repeated this many times evenly around the axis, as in filaments built from several protofilaments"""
        axis: VectorSocket
        """Direction of the helical axis"""
        centre: VectorSocket
        """Point the helical axis passes through"""
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
        count: InputInteger = 10,
        rise: InputFloat = 27.5,
        twist: InputFloat = -2.909464,
        axial_symmetry: InputInteger = 1,
        axis: InputVector = None,
        centre: InputVector = None,
        animate: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Count": count,
                "Rise": rise,
                "Twist": twist,
                "Axial Symmetry": axial_symmetry,
                "Axis": axis,
                "Centre": centre,
                "Animate": animate,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to replicate along the helix"
        )
        count = tree.inputs.integer(
            "Count",
            10,
            description="Number of subunits to place along the helix, each repeated Axial Symmetry times around the axis",
            min_value=1,
            max_value=10000,
        )
        rise = tree.inputs.float(
            "Rise",
            27.5,
            description="Distance advanced along the axis per subunit, in Angstrom",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        twist = tree.inputs.float(
            "Twist",
            -2.909464,
            description="Angle rotated about the axis per subunit",
            subtype="ANGLE",
        )
        axial_symmetry = tree.inputs.integer(
            "Axial Symmetry",
            1,
            description="Cn rotational symmetry about the helical axis: each subunit is repeated this many times evenly around the axis, as in filaments built from several protofilaments",
            min_value=1,
            max_value=100,
        )
        axis = tree.inputs.vector(
            "Axis",
            (0.0, 0.0, 1.0),
            description="Direction of the helical axis",
            subtype="XYZ",
        )
        centre = tree.inputs.vector(
            "Centre",
            (0.0, 0.0, 0.0),
            description="Point the helical axis passes through",
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

        index = g.Index()
        with g.Frame("Helical operators"):
            integer_math = index.o.index // axial_symmetry
            axis_angle_to_rotation = g.AxisAngleToRotation(
                axis=axis,
                angle=twist * integer_math
                + index.o.index % axial_symmetry * (math.tau / axial_symmetry),
            )
            _string = g.String(
                string="Copy k is subunit k div n of the helix, where n is Axial Symmetry: it is twisted that many times about the axis and raised that many times along it, then turned a further (k mod n) / n of a full turn about the axis. Rise is given in Angstrom and converted to world units; the defaults are actin's 27.5 A rise and -166.7 degree twist."
            )
            vector_math = axis.normalize() * (
                AngstromToWorld(angstrom=rise).o.world * integer_math
            )
        (
            SymmetryInstance(
                geometry=geometry,
                points=g.Points(count=count * axial_symmetry),
                rotation=axis_angle_to_rotation,
                translation=vector_math,
                centre=centre,
                animate=animate,
            )
            >> instances
        )


ASSET = SymmetryHelical

ASSET_METADATA = {
    "description": "Build a helical assembly: each subunit advances along the axis by the rise and rotates about it by the twist",
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
