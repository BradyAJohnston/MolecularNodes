# Node-group asset "Filter Render Layer" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    CustomGeometryGroup,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputGeometry


class MN_render_layer_current(CustomGeometryGroup):
    _name = ".MN_render_layer_current"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        layer = tree.outputs.integer(
            "Layer",
            description="Render layer of the view layer being evaluated, driven by Molecular Nodes from the view layer. -1 when it isn't a render layer, which keeps everything",
        )

        integer = g.Integer(-1)

        integer >> layer


class FilterRenderLayer(AssetGeometryGroup):
    """
    Keep only the geometry on the render layer of the view layer being evaluated, as set by Set Render Layer. Keeps everything outside render layers

    Parameters
    ----------
    geometry : InputGeometry
        Geometry

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Filter Render Layer"
    _asset_name = "Filter Render Layer"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Keep only the geometry on the render layer of the view layer being evaluated, as set by Set Render Layer. Keeps everything outside render layers"
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

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
    ):
        super().__init__(**{"Geometry": geometry})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        geometry_1 = tree.outputs.geometry("Geometry")

        mn_render_layer_current = MN_render_layer_current()
        named_attribute = g.NamedAttribute.integer("mn_render_layer")
        boolean_math = (mn_render_layer_current >= 0) & g.Compare.integer.not_equal(
            named_attribute.o.attribute, 0
        )
        compare = g.Compare.integer.not_equal(
            named_attribute.o.attribute, mn_render_layer_current.o.layer + 1
        )
        boolean_math_1 = boolean_math & compare
        separate_components = g.SeparateComponents(geometry)
        join_geometry = g.JoinGeometry(
            (
                separate_components >> g.DeleteGeometry.point(selection=boolean_math_1),
                g.DeleteGeometry.point(separate_components.o.curve, boolean_math_1),
                g.DeleteGeometry.point(
                    separate_components.o.point_cloud, boolean_math_1
                ),
                g.DeleteGeometry.instance(
                    separate_components.o.instances, boolean_math_1
                ),
                separate_components.o.grease_pencil,
                separate_components.o.volume,
            )
        )

        join_geometry >> geometry_1


ASSET = FilterRenderLayer

ASSET_METADATA = {
    "description": "Keep only the geometry on the render layer of the view layer being evaluated, as set by Set Render Layer. Keeps everything outside render layers",
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
