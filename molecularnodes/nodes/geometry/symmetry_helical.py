# Node-group asset "Symmetry Helical" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
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
        Number of subunits to place along the helix
    rise : InputFloat
        Distance advanced along the axis per subunit, in Angstrom
    twist : InputFloat
        Angle rotated about the axis per subunit
    axis : InputVector
        Direction of the helical axis
    centre : InputVector
        Point the helical axis passes through
    factor : InputFloat
        0 places every copy back on the original, 1 builds the full symmetry
    stagger : InputFloat
        Delay each copy by its order, so at 1 the last copy only starts moving as the first finishes

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to replicate along the helix
    i.count : IntegerSocket
        Number of subunits to place along the helix
    i.rise : FloatSocket
        Distance advanced along the axis per subunit, in Angstrom
    i.twist : FloatSocket
        Angle rotated about the axis per subunit
    i.axis : VectorSocket
        Direction of the helical axis
    i.centre : VectorSocket
        Point the helical axis passes through
    i.factor : FloatSocket
        0 places every copy back on the original, 1 builds the full symmetry
    i.stagger : FloatSocket
        Delay each copy by its order, so at 1 the last copy only starts moving as the first finishes

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
        """Number of subunits to place along the helix"""
        rise: FloatSocket
        """Distance advanced along the axis per subunit, in Angstrom"""
        twist: FloatSocket
        """Angle rotated about the axis per subunit"""
        axis: VectorSocket
        """Direction of the helical axis"""
        centre: VectorSocket
        """Point the helical axis passes through"""
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
        count: InputInteger = 10,
        rise: InputFloat = 27.5,
        twist: InputFloat = -2.909464,
        axis: InputVector = None,
        centre: InputVector = None,
        factor: InputFloat = 1.0,
        stagger: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Count": count,
                "Rise": rise,
                "Twist": twist,
                "Axis": axis,
                "Centre": centre,
                "Factor": factor,
                "Stagger": stagger,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to replicate along the helix"
        )
        count = tree.inputs.integer(
            "Count",
            10,
            description="Number of subunits to place along the helix",
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
        with g.Frame("Helical operators"):
            axis_angle_to_rotation = g.AxisAngleToRotation(
                axis=axis, angle=twist * index
            )
            _string = g.String(
                string="Subunit k is twisted k times about the axis and raised k times along it. Rise is given in Angstrom and converted to world units; the defaults are actin's 27.5 A rise and -166.7 degree twist."
            )
            vector_math = axis.normalize() * (
                AngstromToWorld(angstrom=rise).o.world * index
            )
        (
            SymmetryInstance(
                geometry=geometry,
                points=g.Points(count=count, radius=0.1),
                rotation=axis_angle_to_rotation,
                translation=vector_math,
                centre=centre,
                factor=factor,
                stagger=stagger,
            )
            >> instances
        )


ASSET = SymmetryHelical

ASSET_METADATA = {
    "description": "Build a helical assembly: each subunit advances along the axis by the rise and rotates about it by the twist",
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
