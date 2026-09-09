# Node-group asset 'Boolean Run Trim' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class BooleanRunTrim(AssetGeometryGroup):
    """
    Boolean Run Trim

    Parameters
    ----------
    boolean : InputBoolean
        The `Boolean` field to check for continuous runs of `True` values
    start : InputInteger
        The first `n` values of a run of `True` values are made to be `False`
    end : InputInteger
        Ther last `n` values of a run of `True` values become `False`
    size : InputInteger
        A run of `True` values becomes `False` if shorter than this minimum length

    Inputs
    ------
    i.boolean : BooleanSocket
        The `Boolean` field to check for continuous runs of `True` values
    i.start : IntegerSocket
        The first `n` values of a run of `True` values are made to be `False`
    i.end : IntegerSocket
        Ther last `n` values of a run of `True` values become `False`
    i.size : IntegerSocket
        A run of `True` values becomes `False` if shorter than this minimum length

    Outputs
    -------
    o.boolean : BooleanSocket
        The `Boolean` field that has been trimmed
    """

    _name = "Boolean Run Trim"
    _asset_name = "Boolean Run Trim"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.boolean_run_trim"}

    class _Inputs(SocketAccessor):
        boolean: BooleanSocket
        """The `Boolean` field to check for continuous runs of `True` values"""
        start: IntegerSocket
        """The first `n` values of a run of `True` values are made to be `False`"""
        end: IntegerSocket
        """Ther last `n` values of a run of `True` values become `False`"""
        size: IntegerSocket
        """A run of `True` values becomes `False` if shorter than this minimum length"""

    class _Outputs(SocketAccessor):
        boolean: BooleanSocket
        """The `Boolean` field that has been trimmed"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        boolean: InputBoolean = True,
        start: InputInteger = 0,
        end: InputInteger = 0,
        size: InputInteger = 10,
    ):
        super().__init__(
            **{"Boolean": boolean, "Start": start, "End": end, "Size": size}
        )

    def _build_group(self, tree):
        boolean = tree.inputs.boolean(
            "Boolean",
            True,
            description="The `Boolean` field to check for continuous runs of `True` values",
            hide_value=True,
        )
        start = tree.inputs.integer(
            "Start",
            0,
            description="The first `n` values of a run of `True` values are made to be `False`",
            min_value=0,
        )
        end = tree.inputs.integer(
            "End",
            0,
            description="Ther last `n` values of a run of `True` values become `False`",
        )
        size = tree.inputs.integer(
            "Size",
            10,
            description="A run of `True` values becomes `False` if shorter than this minimum length",
            min_value=0,
        )
        boolean_1 = tree.outputs.boolean(
            "Boolean", description="The `Boolean` field that has been trimmed"
        )

        accumulate_field = g.AccumulateField.point.integer(
            group_index=g.AccumulateField.point.integer(~boolean).o.trailing
        )
        (
            (
                boolean
                & (accumulate_field.o.leading > start)
                & (accumulate_field.o.total > size)
                & (accumulate_field.o.total - accumulate_field.o.leading > end)
            )
            >> boolean_1
        )


ASSET = BooleanRunTrim

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
