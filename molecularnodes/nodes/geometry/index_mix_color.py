# Node-group asset "Index Mix Color" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputFloat
from .fractionate_float import FractionateFloat


class IndexMixColor(AssetGeometryGroup):
    """
    Index Mix Color

    Parameters
    ----------
    color : InputColor
        The field to interpolate based on the input `Index`
    index : InputFloat
        The floor and ceiling of this Index value is taken and used for sampling, the fraction of this value is then used to mix between the sampled values

    Inputs
    ------
    i.color : ColorSocket
        The field to interpolate based on the input `Index`
    i.index : FloatSocket
        The floor and ceiling of this Index value is taken and used for sampling, the fraction of this value is then used to mix between the sampled values

    Outputs
    -------
    o.color : ColorSocket
        The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`
    o.from_ : IntegerSocket
        The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`
    o.to : IntegerSocket
        The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`
    """

    _name = "Index Mix Color"
    _asset_name = "Index Mix Color"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.index_mix_color"}

    class _Inputs(SocketAccessor):
        color: ColorSocket
        """The field to interpolate based on the input `Index`"""
        index: FloatSocket
        """The floor and ceiling of this Index value is taken and used for sampling, the fraction of this value is then used to mix between the sampled values"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`"""
        from_: IntegerSocket
        """The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`"""
        to: IntegerSocket
        """The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        color: InputColor = None,
        index: InputFloat = 0.0,
    ):
        super().__init__(**{"Color": color, "Index": index})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        color = tree.inputs.color(
            "Color",
            (1.0, 1.0, 1.0, 1.0),
            description="The field to interpolate based on the input `Index`",
        )
        index = tree.inputs.float(
            "Index",
            0.0,
            description="The floor and ceiling of this Index value is taken and used for sampling, the fraction of this value is then used to mix between the sampled values",
        )
        color_1 = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 1.0),
            description="The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`",
        )
        from_ = tree.outputs.integer(
            "From",
            description="The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`",
        )
        to = tree.outputs.integer(
            "To",
            description="The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`",
        )

        group = FractionateFloat(value=index)
        mix = g.Mix(
            factor_float=group.o.fraction,
            a_color=color.point.at(group.o.floor),
            b_color=color.point.at(group.o.ceiling),
            data_type="RGBA",
            clamp_factor=True,
            clamp_result=True,
        )

        mix.o.result_color >> color_1
        group.o.floor >> from_
        group.o.ceiling >> to


ASSET = IndexMixColor

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
