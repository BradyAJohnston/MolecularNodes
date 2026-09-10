# Node-group asset 'Transform Accumulate' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    MatrixSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputInteger, InputMatrix, InputMenu
from ._shared.accumulate_domain_transform import AccumulateDomainTransform


class TransformAccumulate(AssetGeometryGroup):
    """
    Transform Accumulate

    Parameters
    ----------
    domain : InputMenu | Literal["Point", "Edge", "Face", "Face Corner", "Spline", "Instance"]
        Domain on which to accumulate transforms
    accumulate : InputBoolean
        Selected transforms are included in the accumulation, non-selected transforms become identity matrices
    transform : InputMatrix
        Transform field to accumulate
    group_id : InputInteger
        Transform field is accumulated individually for each `Group ID`

    Inputs
    ------
    i.domain : MenuSocket
        Domain on which to accumulate transforms
    i.accumulate : BooleanSocket
        Selected transforms are included in the accumulation, non-selected transforms become identity matrices
    i.transform : MatrixSocket
        Transform field to accumulate
    i.group_id : IntegerSocket
        Transform field is accumulated individually for each `Group ID`

    Outputs
    -------
    o.transform : MatrixSocket
        The accumulating Transform field
    """

    _name = "Transform Accumulate"
    _asset_name = "Transform Accumulate"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.transform_accumulate"}

    class _Inputs(SocketAccessor):
        domain: MenuSocket
        """Domain on which to accumulate transforms"""
        accumulate: BooleanSocket
        """Selected transforms are included in the accumulation, non-selected transforms become identity matrices"""
        transform: MatrixSocket
        """Transform field to accumulate"""
        group_id: IntegerSocket
        """Transform field is accumulated individually for each `Group ID`"""

    class _Outputs(SocketAccessor):
        transform: MatrixSocket
        """The accumulating Transform field"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        domain: InputMenu
        | Literal[
            "Point", "Edge", "Face", "Face Corner", "Spline", "Instance"
        ] = "Point",
        accumulate: InputBoolean = True,
        transform: InputMatrix = None,
        group_id: InputInteger = 0,
    ):
        super().__init__(
            **{
                "Domain": domain,
                "Accumulate": accumulate,
                "Transform": transform,
                "Group ID": group_id,
            }
        )

    def _build_group(self, tree):
        domain = tree.inputs.menu(
            "Domain",
            description="Domain on which to accumulate transforms",
            optional_label=True,
        )
        accumulate = tree.inputs.boolean(
            "Accumulate",
            True,
            description="Selected transforms are included in the accumulation, non-selected transforms become identity matrices",
            hide_value=True,
        )
        transform = tree.inputs.matrix(
            "Transform", description="Transform field to accumulate", hide_value=True
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="Transform field is accumulated individually for each `Group ID`",
            hide_value=True,
        )
        transform_1 = tree.outputs.matrix(
            "Transform", description="The accumulating Transform field"
        )

        (
            AccumulateDomainTransform(
                domain=domain,
                transform=accumulate.switch.matrix(true=transform),
                group_id=group_id,
            )
            >> transform_1
        )

        domain.default_value = "Point"


ASSET = TransformAccumulate

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
