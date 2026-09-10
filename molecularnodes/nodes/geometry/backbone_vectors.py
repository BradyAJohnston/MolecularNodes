# Node-group asset "Backbone Vectors" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputMenu
from .backbone_positions import BackbonePositions


class BackboneVectors(AssetGeometryGroup):
    """
    The `Vectors` that are useful for a curve when reading from a peptide backbone

    Parameters
    ----------
    method : InputMenu | Literal["Read", "Compute"]
        Method

    Inputs
    ------
    i.method : MenuSocket
        Method

    Outputs
    -------
    o.normal : VectorSocket
        The vector used for the `Normal` of a curve when reading positions from a peptide backbone
    o.tangent : VectorSocket
        The vector used as the `Tangent` for a curve when reading values from a peptide backbone
    o.bitangent : VectorSocket
        Cross product of the Normal and Tangent
    """

    _name = "Backbone Vectors"
    _asset_name = "Backbone Vectors"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "The `Vectors` that are useful for a curve when reading from a peptide backbone",
        "node_tool_idname": "geometry.backbone_vectors",
    }

    class _Inputs(SocketAccessor):
        method: MenuSocket
        """Method"""

    class _Outputs(SocketAccessor):
        normal: VectorSocket
        """The vector used for the `Normal` of a curve when reading positions from a peptide backbone"""
        tangent: VectorSocket
        """The vector used as the `Tangent` for a curve when reading values from a peptide backbone"""
        bitangent: VectorSocket
        """Cross product of the Normal and Tangent"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        method: InputMenu | Literal["Read", "Compute"] = "Compute",
    ):
        super().__init__(**{"Method": method})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        method = tree.inputs.menu("Method", expanded=True, optional_label=True)
        normal = tree.outputs.vector(
            "Normal",
            description="The vector used for the `Normal` of a curve when reading positions from a peptide backbone",
        )
        tangent = tree.outputs.vector(
            "Tangent",
            description="The vector used as the `Tangent` for a curve when reading values from a peptide backbone",
        )
        bitangent = tree.outputs.vector(
            "Bitangent", description="Cross product of the Normal and Tangent"
        )

        group = BackbonePositions(method=method)
        mix = g.Mix(
            a_vector=group.o.c,
            b_vector=group.o.n,
            factor_float=0.45,
            data_type="VECTOR",
            clamp_factor=True,
        )
        vector_math = (group.o.c - group.o.n).normalize()
        vector_math_1 = (mix.o.result_vector - group.o.ca).normalize()
        vector_math.cross(vector_math_1).normalize() >> bitangent

        vector_math_1 >> normal
        vector_math >> tangent

        method.default_value = "Compute"


ASSET = BackboneVectors

ASSET_METADATA = {
    "description": "The `Vectors` that are useful for a curve when reading from a peptide backbone",
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
