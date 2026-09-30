# Node-group asset "Fade Geometry" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputGeometry


class FadeGeometry(AssetGeometryGroup):
    """
    Fade Geometry

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to fade, such as the output of a style node
    fade : InputFloat
        Multiplier for the alpha of the `Color` attribute. At `1` the geometry is passed through unchanged, at `0` it is removed

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to fade, such as the output of a style node
    i.fade : FloatSocket
        Multiplier for the alpha of the `Color` attribute. At `1` the geometry is passed through unchanged, at `0` it is removed

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry with the faded alpha stored in `Color`
    """

    _name = "Fade Geometry"
    _asset_name = "Fade Geometry"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to fade, such as the output of a style node"""
        fade: FloatSocket
        """Multiplier for the alpha of the `Color` attribute. At `1` the geometry is passed through unchanged, at `0` it is removed"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry with the faded alpha stored in `Color`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        fade: InputFloat = 1.0,
    ):
        super().__init__(**{"Geometry": geometry, "Fade": fade})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry",
            description="Geometry to fade, such as the output of a style node",
        )
        fade = tree.inputs.float(
            "Fade",
            1.0,
            description="Multiplier for the alpha of the `Color` attribute. At `1` the geometry is passed through unchanged, at `0` it is removed",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        geometry_1 = tree.outputs.geometry(
            "Geometry", description="Geometry with the faded alpha stored in `Color`"
        )

        get_geometry_component = g.GetGeometryComponent(
            geometry=geometry, type="Instances"
        )
        attribute = g.NamedAttribute.color("Color").o.attribute
        combine_color = g.CombineColor(
            red=attribute.r,
            green=attribute.g,
            blue=attribute.b,
            alpha=attribute.a * fade,
        )
        store_named_attribute = (
            get_geometry_component
            >> g.StoreNamedAttribute.point.color(name="Color", value=combine_color)
        )
        store_named_attribute_1 = g.StoreNamedAttribute.instance.color(
            get_geometry_component.o.component, name="Color", value=combine_color
        )
        switch = (fade > 0.0).switch.geometry(
            true=g.JoinGeometry(
                geometry=(store_named_attribute, store_named_attribute_1)
            )
        )
        (fade < 1.0).switch.geometry(geometry, switch) >> geometry_1


ASSET = FadeGeometry

ASSET_METADATA = {
    "description": "Fade geometry in and out by multiplying the alpha of its `Color` attribute on the point and instance domains. A `Fade` of `0` or below returns empty geometry, and a `Fade` of `1` or above passes the geometry through without any computation.",
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
