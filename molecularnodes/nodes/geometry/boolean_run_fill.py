# Node-group asset 'Boolean Run Fill' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputInteger


class BooleanRunFill(AssetGeometryGroup):
    """
    Moving down the points, fill in `False` values that are equal to or less than the `Fill Size`

    Parameters
    ----------
    boolean : InputBoolean
        The `Boolean` field to potentially fill gaps of `False` with
    fill_size : InputInteger
        A run of `False` values equal to or less than this size will become `True`

    Inputs
    ------
    i.boolean : BooleanSocket
        The `Boolean` field to potentially fill gaps of `False` with
    i.fill_size : IntegerSocket
        A run of `False` values equal to or less than this size will become `True`

    Outputs
    -------
    o.boolean : BooleanSocket
        The `Boolean` array with gaps potentially filled
    """

    _name = "Boolean Run Fill"
    _asset_name = "Boolean Run Fill"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Moving down the points, fill in `False` values that are equal to or less than the `Fill Size`",
        "node_tool_idname": "geometry.boolean_run_fill",
    }

    class _Inputs(SocketAccessor):
        boolean: BooleanSocket
        """The `Boolean` field to potentially fill gaps of `False` with"""
        fill_size: IntegerSocket
        """A run of `False` values equal to or less than this size will become `True`"""

    class _Outputs(SocketAccessor):
        boolean: BooleanSocket
        """The `Boolean` array with gaps potentially filled"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        boolean: InputBoolean = True,
        fill_size: InputInteger = 3,
    ):
        super().__init__(**{"Boolean": boolean, "Fill Size": fill_size})

    def _build_group(self, tree):
        boolean = tree.inputs.boolean(
            "Boolean",
            True,
            description="The `Boolean` field to potentially fill gaps of `False` with",
            hide_value=True,
        )
        fill_size = tree.inputs.integer(
            "Fill Size",
            3,
            description="A run of `False` values equal to or less than this size will become `True`",
        )
        boolean_1 = tree.outputs.boolean(
            "Boolean", description="The `Boolean` array with gaps potentially filled"
        )

        accumulate_field = g.AccumulateField.point.integer(
            group_index=g.AccumulateField.point.integer(boolean).o.leading
        )
        (
            (
                boolean
                | (accumulate_field.o.trailing <= fill_size)
                & (accumulate_field.o.total <= fill_size)
            )
            >> boolean_1
        )


ASSET = BooleanRunFill

ASSET_METADATA = {
    "description": "Moving down the points, fill in `False` values that are equal to or less than the `Fill Size`",
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
