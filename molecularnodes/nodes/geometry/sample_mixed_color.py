# Node-group asset "Sample Mixed Color" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputFloat, InputGeometry
from .index_mix_color import IndexMixColor


class SampleMixedColor(AssetGeometryGroup):
    """
    Sample Mixed Color

    Parameters
    ----------
    geometry : InputGeometry
        The geometry to sample the values from
    color : InputColor
        The field to mix and evaluate on the sample geometry
    index : InputFloat
        The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes

    Inputs
    ------
    i.geometry : GeometrySocket
        The geometry to sample the values from
    i.color : ColorSocket
        The field to mix and evaluate on the sample geometry
    i.index : FloatSocket
        The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes

    Outputs
    -------
    o.color : ColorSocket
        The evaluated and mixed field, sampled from the sample geometry at the given `Index`
    """

    _name = "Sample Mixed Color"
    _asset_name = "Sample Mixed Color"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.sample_mixed_color"}

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """The geometry to sample the values from"""
        color: ColorSocket
        """The field to mix and evaluate on the sample geometry"""
        index: FloatSocket
        """The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The evaluated and mixed field, sampled from the sample geometry at the given `Index`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        color: InputColor = None,
        index: InputFloat = 0.0,
    ):
        super().__init__(**{"Geometry": geometry, "Color": color, "Index": index})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="The geometry to sample the values from"
        )
        color = tree.inputs.color(
            "Color",
            (0.0, 0.0, 0.0, 1.0),
            description="The field to mix and evaluate on the sample geometry",
            hide_value=True,
        )
        index = tree.inputs.float(
            "Index",
            0.0,
            description="The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes",
        )
        color_1 = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 1.0),
            description="The evaluated and mixed field, sampled from the sample geometry at the given `Index`",
        )

        group = IndexMixColor(color=color, index=index)
        (
            geometry
            >> g.SampleIndex(
                value=group.o.color,
                index=group.o.from_,
                data_type="FLOAT_COLOR",
                clamp=True,
            )
            >> color_1
        )


ASSET = SampleMixedColor

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
