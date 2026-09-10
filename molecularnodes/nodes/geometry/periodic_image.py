# Node-group asset "Periodic Image" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class PeriodicImage(AssetGeometryGroup):
    """
    Periodic Image

    Parameters
    ----------
    a : InputVector
        A
    b : InputVector
        B
    c : InputVector
        C
    image_a : InputInteger
        Image A
    image_b : InputInteger
        Image B
    image_c : InputInteger
        Image C

    Inputs
    ------
    i.a : VectorSocket
        A
    i.b : VectorSocket
        B
    i.c : VectorSocket
        C
    i.image_a : IntegerSocket
        Image A
    i.image_b : IntegerSocket
        Image B
    i.image_c : IntegerSocket
        Image C

    Outputs
    -------
    o.points : GeometrySocket
        Points
    """

    _name = "Periodic Image"
    _asset_name = "Periodic Image"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        a: VectorSocket
        """A"""
        b: VectorSocket
        """B"""
        c: VectorSocket
        """C"""
        image_a: IntegerSocket
        """Image A"""
        image_b: IntegerSocket
        """Image B"""
        image_c: IntegerSocket
        """Image C"""

    class _Outputs(SocketAccessor):
        points: GeometrySocket
        """Points"""

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
        image_a: InputInteger = 0,
        image_b: InputInteger = 0,
        image_c: InputInteger = 0,
    ):
        super().__init__(
            **{
                "A": a,
                "B": b,
                "C": c,
                "Image A": image_a,
                "Image B": image_b,
                "Image C": image_c,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        a = tree.inputs.vector(
            "A", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0, subtype="XYZ"
        )
        b = tree.inputs.vector(
            "B", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0, subtype="XYZ"
        )
        c_ = tree.inputs.vector(
            "C", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0, subtype="XYZ"
        )
        with tree.inputs.panel("Image"):
            image_a = tree.inputs.integer("Image A", 0)
            image_b = tree.inputs.integer("Image B", 0)
            image_c = tree.inputs.integer("Image C", 0)
        points = tree.outputs.geometry("Points")

        points_1 = g.Points(
            position=a * image_a + b * image_b + c_ * image_c, radius=0.1
        )

        points_1 >> points


ASSET = PeriodicImage

ASSET_METADATA = {
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
