# Node-group asset 'Primitive Gimbal' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
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
from .primitive_arrow import PrimitiveArrow


class PrimitiveGimbal(AssetGeometryGroup):
    """
    Primitive Gimbal

    Parameters
    ----------
    vertices : InputInteger
        Vertices
    x : InputColor
        X
    y : InputColor
        Y
    z : InputColor
        Z
    material : InputMaterial
        Material to apply to the resulting geometry
    height : InputFloat
        Height

    Inputs
    ------
    i.vertices : IntegerSocket
        Vertices
    i.x : ColorSocket
        X
    i.y : ColorSocket
        Y
    i.z : ColorSocket
        Z
    i.material : MaterialSocket
        Material to apply to the resulting geometry
    i.height : FloatSocket
        Height

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Primitive Gimbal"
    _asset_name = "Primitive Gimbal"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.primitive_gimbal"}

    class _Inputs(SocketAccessor):
        vertices: IntegerSocket
        """Vertices"""
        x: ColorSocket
        """X"""
        y: ColorSocket
        """Y"""
        z: ColorSocket
        """Z"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""
        height: FloatSocket
        """Height"""

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
        x: InputColor = None,
        y: InputColor = None,
        z: InputColor = None,
        material: InputMaterial = None,
        height: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Vertices": vertices,
                "X": x,
                "Y": y,
                "Z": z,
                "Material": material,
                "Height": height,
            }
        )

    def _build_group(self, tree):
        vertices = tree.inputs.integer("Vertices", 6, min_value=3, max_value=16)
        x = tree.inputs.color("X", (0.623968, 0.01033, 0.063011, 1.0))
        y = tree.inputs.color("Y", (0.076185, 0.623968, 0.084375, 1.0))
        z = tree.inputs.color("Z", (0.0, 0.00091, 0.623968, 1.0))
        material = tree.inputs.material(
            "Material", description="Material to apply to the resulting geometry"
        )
        height = tree.inputs.float(
            "Height", 1.0, min_value=0.0, max_value=10_000.0, subtype="DISTANCE"
        )
        geometry = tree.outputs.geometry("Geometry")

        repeat_zone = g.RepeatZone(3)
        geometry_1 = repeat_zone.items.geometry("Geometry")
        integer = repeat_zone.items.integer("Integer")
        index_switch = g.IndexSwitch.rotation(
            integer.current,
            ((0.0, math.pi / 2, 0.0), (-math.pi / 2, 0.0, 0.0), (0.0, 0.0, 0.0)),
        )
        transform_geometry = PrimitiveArrow(
            vertices=vertices,
            height=height,
            value=g.IndexSwitch.color(integer.current, (x, y, z)),
            material=material,
        ) >> g.TransformGeometry(rotation=index_switch)
        (
            g.JoinGeometry(geometry=(geometry_1.current, transform_geometry))
            >> geometry_1.next
        )
        integer.current + 1.0 >> integer.next

        geometry_1.result >> geometry


ASSET = PrimitiveGimbal

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
