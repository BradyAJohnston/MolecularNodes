# Node-group asset "Lattice Grid" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputInteger, InputVector


class LatticeGrid(AssetGeometryGroup):
    """
    Lattice Grid

    Parameters
    ----------
    a : InputVector
        Box vector A in world units, such as from Periodic Box
    b : InputVector
        Box vector B in world units, such as from Periodic Box
    c : InputVector
        Box vector C in world units, such as from Periodic Box
    x : InputInteger
        Number of lattice points along box vector A
    y : InputInteger
        Number of lattice points along box vector B
    z : InputInteger
        Number of lattice points along box vector C

    Inputs
    ------
    i.a : VectorSocket
        Box vector A in world units, such as from Periodic Box
    i.b : VectorSocket
        Box vector B in world units, such as from Periodic Box
    i.c : VectorSocket
        Box vector C in world units, such as from Periodic Box
    i.x : IntegerSocket
        Number of lattice points along box vector A
    i.y : IntegerSocket
        Number of lattice points along box vector B
    i.z : IntegerSocket
        Number of lattice points along box vector C

    Outputs
    -------
    o.grid : GeometrySocket
        A point at each lattice point, starting at the origin and going along box vectors A, B and C
    """

    _name = "Lattice Grid"
    _asset_name = "Lattice Grid"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        a: VectorSocket
        """Box vector A in world units, such as from Periodic Box"""
        b: VectorSocket
        """Box vector B in world units, such as from Periodic Box"""
        c: VectorSocket
        """Box vector C in world units, such as from Periodic Box"""
        x: IntegerSocket
        """Number of lattice points along box vector A"""
        y: IntegerSocket
        """Number of lattice points along box vector B"""
        z: IntegerSocket
        """Number of lattice points along box vector C"""

    class _Outputs(SocketAccessor):
        grid: GeometrySocket
        """A point at each lattice point, starting at the origin and going along box vectors A, B and C"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        a: InputVector = None,
        b: InputVector = None,
        c: InputVector = None,
        x: InputInteger = 3,
        y: InputInteger = 3,
        z: InputInteger = 3,
    ):
        super().__init__(**{"A": a, "B": b, "C": c, "X": x, "Y": y, "Z": z})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        a = tree.inputs.vector(
            "A",
            (0.0, 0.0, 1.0),
            description="Box vector A in world units, such as from Periodic Box",
            subtype="XYZ",
        )
        b = tree.inputs.vector(
            "B",
            (0.0, 0.0, 1.0),
            description="Box vector B in world units, such as from Periodic Box",
            subtype="XYZ",
        )
        c_ = tree.inputs.vector(
            "C",
            (0.0, 0.0, 1.0),
            description="Box vector C in world units, such as from Periodic Box",
            subtype="XYZ",
        )
        with tree.inputs.panel("Count"):
            x = tree.inputs.integer(
                "X",
                3,
                description="Number of lattice points along box vector A",
                min_value=1,
                max_value=10000,
            )
            y = tree.inputs.integer(
                "Y",
                3,
                description="Number of lattice points along box vector B",
                min_value=1,
                max_value=10000,
            )
            z = tree.inputs.integer(
                "Z",
                3,
                description="Number of lattice points along box vector C",
                min_value=1,
                max_value=10000,
            )
        grid = tree.outputs.geometry(
            "Grid",
            description="A point at each lattice point, starting at the origin and going along box vectors A, B and C",
        )

        (
            g.MeshLine(count=x, offset=a)
            >> g.InstanceOnPoints(instance=g.MeshLine(count=y, offset=b))
            >> g.InstanceOnPoints(instance=g.MeshLine(count=z, offset=c_))
            >> g.RealizeInstances(realize_to_point_domain=True)
            >> grid
        )


ASSET = LatticeGrid

ASSET_METADATA = {
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
