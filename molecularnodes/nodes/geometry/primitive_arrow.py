# Node-group asset "Primitive Arrow" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MaterialSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputFloat, InputInteger, InputMaterial


class PrimitiveArrow(AssetGeometryGroup):
    """
    Primitive Arrow

    Parameters
    ----------
    vertices : InputInteger
        Vertices
    height : InputFloat
        Height
    ratio : InputFloat
        Ratio
    value : InputColor
        Value
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.vertices : IntegerSocket
        Vertices
    i.height : FloatSocket
        Height
    i.ratio : FloatSocket
        Ratio
    i.value : ColorSocket
        Value
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Primitive Arrow"
    _asset_name = "Primitive Arrow"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.primitive_arrow"}

    class _Inputs(SocketAccessor):
        vertices: IntegerSocket
        """Vertices"""
        height: FloatSocket
        """Height"""
        ratio: FloatSocket
        """Ratio"""
        value: ColorSocket
        """Value"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        vertices: InputInteger = 6,
        height: InputFloat = 1.0,
        ratio: InputFloat = 0.7,
        value: InputColor = None,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Vertices": vertices,
                "Height": height,
                "Ratio": ratio,
                "Value": value,
                "Material": material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        vertices = tree.inputs.integer("Vertices", 6, min_value=3, max_value=16)
        height = tree.inputs.float(
            "Height", 1.0, min_value=0.0, max_value=10_000.0, subtype="DISTANCE"
        )
        ratio = tree.inputs.float(
            "Ratio", 0.7, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        value = tree.inputs.color("Value", (0.8, 0.8, 0.8, 1.0))
        material = tree.inputs.material(
            "Material", description="Material to apply to the resulting geometry"
        )
        geometry = tree.outputs.geometry("Geometry")

        math_1 = height * ratio
        math_2 = height - math_1
        transform_geometry = g.TransformGeometry(
            geometry=g.Cone(vertices=vertices, depth=math_2, radius_bottom=0.2),
            translation=g.CombineXYZ(z=height - math_2),
        )
        transform_geometry_1 = g.TransformGeometry(
            geometry=g.Cylinder(vertices=vertices, depth=math_1, radius=0.08),
            translation=g.CombineXYZ(z=math_1 / 2.0),
        )
        (
            g.JoinGeometry(geometry=(transform_geometry, transform_geometry_1))
            >> g.StoreNamedAttribute.point.color(name="Color", value=value)
            >> g.SetMaterial(material=material)
            >> geometry
        )


ASSET = PrimitiveArrow

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
