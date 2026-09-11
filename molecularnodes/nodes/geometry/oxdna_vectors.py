# Node-group asset "oxDNA Vectors" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputInteger
from .angstrom_to_world import AngstromToWorld


class OxDNAVectors(AssetGeometryGroup):
    """
    oxDNA Vectors

    Parameters
    ----------
    index : InputInteger
        Index

    Inputs
    ------
    i.index : IntegerSocket
        Index

    Outputs
    -------
    o.base_vector : VectorSocket
        base_vector
    o.base_normal : VectorSocket
        base_normal
    o.backbone_offset : VectorSocket
        Backbone Offset
    o.base_offset : VectorSocket
        Base Offset
    o.rotation : RotationSocket
        Rotation
    """

    _name = "oxDNA Vectors"
    _asset_name = "oxDNA Vectors"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        base_vector: VectorSocket
        base_normal: VectorSocket
        backbone_offset: VectorSocket
        """Backbone Offset"""
        base_offset: VectorSocket
        """Base Offset"""
        rotation: RotationSocket
        """Rotation"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
    ):
        super().__init__(**{"Index": index})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        index = tree.inputs.integer(
            "Index", 0, min_value=0, hide_value=True, default_input="INDEX"
        )
        base_vector = tree.outputs.vector("base_vector")
        base_normal = tree.outputs.vector("base_normal")
        backbone_offset = tree.outputs.vector("Backbone Offset")
        base_offset = tree.outputs.vector("Base Offset")
        rotation = tree.outputs.rotation("Rotation")

        group = AngstromToWorld(angstrom=3.4)
        evaluate_at_index = g.NamedAttribute.vector("base_vector").o.attribute.point.at(
            index
        )
        evaluate_at_index_1 = g.NamedAttribute.vector(
            "base_normal"
        ).o.attribute.point.at(index)
        axes_to_rotation = g.AxesToRotation(
            primary_axis=evaluate_at_index, secondary_axis=evaluate_at_index_1
        )
        evaluate_at_index * group >> base_offset
        (
            evaluate_at_index * (group.o.world * -1.0)
            + evaluate_at_index.cross(evaluate_at_index_1) * group
            >> backbone_offset
        )

        evaluate_at_index >> base_vector
        evaluate_at_index_1 >> base_normal
        axes_to_rotation >> rotation


ASSET = OxDNAVectors

ASSET_METADATA = {
    "description": "Vectors and rotation relevant to the oxDNA model of DNA",
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}
