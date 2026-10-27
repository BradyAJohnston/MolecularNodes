# Node-group asset "Set Render Layer" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputGeometry, InputInteger


class SetRenderLayer(AssetGeometryGroup):
    """
    Put the selected geometry on a Molecular Nodes render layer, stored in the mn_render_layer attribute as layer + 1. Geometry that is never set (0) renders on every layer

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    selection : InputBoolean
        Points and instances to put on the layer. Unselected elements keep their current layer
    layer : InputInteger
        Render layer to put the selection on. Higher layers are composited over lower ones

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.selection : BooleanSocket
        Points and instances to put on the layer. Unselected elements keep their current layer
    i.layer : IntegerSocket
        Render layer to put the selection on. Higher layers are composited over lower ones

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Set Render Layer"
    _asset_name = "Set Render Layer"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "ATTRIBUTE"
    _tree_properties = {
        "description": "Put the selected geometry on a Molecular Nodes render layer, stored in the mn_render_layer attribute as layer + 1. Geometry that is never set (0) renders on every layer"
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        selection: BooleanSocket
        """Points and instances to put on the layer. Unselected elements keep their current layer"""
        layer: IntegerSocket
        """Render layer to put the selection on. Higher layers are composited over lower ones"""

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
        selection: InputBoolean = True,
        layer: InputInteger = 0,
    ):
        super().__init__(
            **{"Geometry": geometry, "Selection": selection, "Layer": layer}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Points and instances to put on the layer. Unselected elements keep their current layer",
            hide_value=True,
        )
        layer = tree.inputs.integer(
            "Layer",
            0,
            description="Render layer to put the selection on. Higher layers are composited over lower ones",
            min_value=0,
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        integer_math = layer + 1
        (
            geometry
            >> g.StoreNamedAttribute.point.integer(
                selection=selection, name="mn_render_layer", value=integer_math
            )
            >> g.StoreNamedAttribute.instance.integer(
                selection=selection, name="mn_render_layer", value=integer_math
            )
            >> geometry_1
        )


ASSET = SetRenderLayer

ASSET_METADATA = {
    "description": "Put the selected geometry on a Molecular Nodes render layer, stored in the mn_render_layer attribute as layer + 1. Geometry that is never set (0) renders on every layer",
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
