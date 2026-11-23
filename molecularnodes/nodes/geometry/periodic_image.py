# Node-group asset "Periodic Image" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class PeriodicImage(AssetGeometryGroup):
    """
    Periodic Image

    Parameters
    ----------
    a : InputVector
        Box vector A in world units, such as from Periodic Box
    b : InputVector
        Box vector B in world units, such as from Periodic Box
    c : InputVector
        Box vector C in world units, such as from Periodic Box
    image_a : InputInteger
        Number of box vectors A from the origin to the image
    image_b : InputInteger
        Number of box vectors B from the origin to the image
    image_c : InputInteger
        Number of box vectors C from the origin to the image

    Inputs
    ------
    i.a : VectorSocket
        Box vector A in world units, such as from Periodic Box
    i.b : VectorSocket
        Box vector B in world units, such as from Periodic Box
    i.c : VectorSocket
        Box vector C in world units, such as from Periodic Box
    i.image_a : IntegerSocket
        Number of box vectors A from the origin to the image
    i.image_b : IntegerSocket
        Number of box vectors B from the origin to the image
    i.image_c : IntegerSocket
        Number of box vectors C from the origin to the image

    Outputs
    -------
    o.points : GeometrySocket
        A single point at the image, Image A × A + Image B × B + Image C × C
    """

    _name = "Periodic Image"
    _asset_name = "Periodic Image"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        a: VectorSocket
        """Box vector A in world units, such as from Periodic Box"""
        b: VectorSocket
        """Box vector B in world units, such as from Periodic Box"""
        c: VectorSocket
        """Box vector C in world units, such as from Periodic Box"""
        image_a: IntegerSocket
        """Number of box vectors A from the origin to the image"""
        image_b: IntegerSocket
        """Number of box vectors B from the origin to the image"""
        image_c: IntegerSocket
        """Number of box vectors C from the origin to the image"""

    class _Outputs(SocketAccessor):
        points: GeometrySocket
        """A single point at the image, Image A × A + Image B × B + Image C × C"""

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
            "A",
            (0.0, 0.0, 0.0),
            description="Box vector A in world units, such as from Periodic Box",
            min_value=-10_000.0,
            max_value=10_000.0,
            subtype="XYZ",
        )
        b = tree.inputs.vector(
            "B",
            (0.0, 0.0, 0.0),
            description="Box vector B in world units, such as from Periodic Box",
            min_value=-10_000.0,
            max_value=10_000.0,
            subtype="XYZ",
        )
        c_ = tree.inputs.vector(
            "C",
            (0.0, 0.0, 0.0),
            description="Box vector C in world units, such as from Periodic Box",
            min_value=-10_000.0,
            max_value=10_000.0,
            subtype="XYZ",
        )
        with tree.inputs.panel("Image"):
            image_a = tree.inputs.integer(
                "Image A",
                0,
                description="Number of box vectors A from the origin to the image",
            )
            image_b = tree.inputs.integer(
                "Image B",
                0,
                description="Number of box vectors B from the origin to the image",
            )
            image_c = tree.inputs.integer(
                "Image C",
                0,
                description="Number of box vectors C from the origin to the image",
            )
        points = tree.outputs.geometry(
            "Points",
            description="A single point at the image, Image A × A + Image B × B + Image C × C",
        )

        points_1 = g.Points(position=a * image_a + b * image_b + c_ * image_c)

        points_1 >> points


ASSET = PeriodicImage

ASSET_METADATA = {
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
