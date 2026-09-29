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
        The generated geometry for the style node group
    fade : InputFloat
        Fade

    Inputs
    ------
    i.geometry : GeometrySocket
        The generated geometry for the style node group
    i.fade : FloatSocket
        Fade

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Fade Geometry"
    _asset_name = "Fade Geometry"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """The generated geometry for the style node group"""
        fade: FloatSocket
        """Fade"""

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
        geometry: InputGeometry = None,
        fade: InputFloat = 1.0,
    ):
        super().__init__(**{"Geometry": geometry, "Fade": fade})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="The generated geometry for the style node group"
        )
        fade = tree.inputs.float(
            "Fade", 1.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        geometry_1 = tree.outputs.geometry("Geometry")

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
    "description": "Set the alpha of the `Color` attribute on the point and instance domains of the geometry. If fade is `0` or below then an empty geometry is returned. If `Fade` is `1.0` or above then computation is skipped.",
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
