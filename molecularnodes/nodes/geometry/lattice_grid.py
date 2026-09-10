# Node-group asset "Lattice Grid" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
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
        A
    b : InputVector
        B
    c : InputVector
        C
    x : InputInteger
        Number of images along the a vector
    y : InputInteger
        Number of images along the b vector
    z : InputInteger
        Number of images along the c vector

    Inputs
    ------
    i.a : VectorSocket
        A
    i.b : VectorSocket
        B
    i.c : VectorSocket
        C
    i.x : IntegerSocket
        Number of images along the a vector
    i.y : IntegerSocket
        Number of images along the b vector
    i.z : IntegerSocket
        Number of images along the c vector

    Outputs
    -------
    o.grid : GeometrySocket
        Grid
    """

    _name = "Lattice Grid"
    _asset_name = "Lattice Grid"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        a: VectorSocket
        """A"""
        b: VectorSocket
        """B"""
        c: VectorSocket
        """C"""
        x: IntegerSocket
        """Number of images along the a vector"""
        y: IntegerSocket
        """Number of images along the b vector"""
        z: IntegerSocket
        """Number of images along the c vector"""

    class _Outputs(SocketAccessor):
        grid: GeometrySocket
        """Grid"""

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
        a = tree.inputs.vector("A", (0.0, 0.0, 1.0), subtype="XYZ")
        b = tree.inputs.vector("B", (0.0, 0.0, 1.0), subtype="XYZ")
        c_ = tree.inputs.vector("C", (0.0, 0.0, 1.0), subtype="XYZ")
        with tree.inputs.panel("Count"):
            x = tree.inputs.integer(
                "X",
                3,
                description="Number of images along the a vector",
                min_value=1,
                max_value=10000,
            )
            y = tree.inputs.integer(
                "Y",
                3,
                description="Number of images along the b vector",
                min_value=1,
                max_value=10000,
            )
            z = tree.inputs.integer(
                "Z",
                3,
                description="Number of images along the c vector",
                min_value=1,
                max_value=10000,
            )
        grid = tree.outputs.geometry("Grid")

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
