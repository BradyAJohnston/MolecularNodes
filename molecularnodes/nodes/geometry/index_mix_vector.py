# Node-group asset "Index Mix Vector" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputVector
from .fractionate_float import FractionateFloat


class IndexMixVector(AssetGeometryGroup):
    """
    Index Mix Vector

    Parameters
    ----------
    value : InputVector
        The field to interpolate based on the input `Index`
    index : InputFloat
        The floor and ceiling of this Index value is taken and used for sampling, the fraction of this value is then used to mix between the sampled values

    Inputs
    ------
    i.value : VectorSocket
        The field to interpolate based on the input `Index`
    i.index : FloatSocket
        The floor and ceiling of this Index value is taken and used for sampling, the fraction of this value is then used to mix between the sampled values

    Outputs
    -------
    o.value : VectorSocket
        The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`
    o.from_ : IntegerSocket
        The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`
    o.to : IntegerSocket
        The mixed value of the field, first evaluating the field at the `From` and `To` Indices then mixing between them based on the fraction of the input `Index`
    """

    _name = "Index Mix Vector"
    _asset_name = "Index Mix Vector"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.index_mix_vector"}

    class _Inputs(SocketAccessor):
        value: VectorSocket
        """The field to interpolate based on the input `Index`"""
        index: FloatSocket
        """The floor and ceiling of this Index value is taken and used for sampling, the fraction of this value is then used to mix between the sampled values"""

    class _Outputs(SocketAccessor):
        value: VectorSocket
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
        value: InputVector = None,
        index: InputFloat = 0.0,
    ):
        super().__init__(**{"Value": value, "Index": index})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        value = tree.inputs.vector(
            "Value",
            (0.0, 0.0, 0.0),
            description="The field to interpolate based on the input `Index`",
        )
        index = tree.inputs.float(
            "Index",
            0.0,
            description="The floor and ceiling of this Index value is taken and used for sampling, the fraction of this value is then used to mix between the sampled values",
        )
        value_1 = tree.outputs.vector(
            "Value",
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
            a_vector=value.point.at(group.o.floor),
            b_vector=value.point.at(group.o.ceiling),
            data_type="VECTOR",
            clamp_factor=True,
        )

        mix.o.result_vector >> value_1
        group.o.floor >> from_
        group.o.ceiling >> to


ASSET = IndexMixVector

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
