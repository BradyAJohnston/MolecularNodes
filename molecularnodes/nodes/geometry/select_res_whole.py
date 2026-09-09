# Node-group asset 'Select Res Whole' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean
from .boolean_any import BooleanAny
from .unique_residue_id import UniqueResidueID


class SelectResWhole(AssetGeometryGroup):
    """
    Select Res Whole

    Parameters
    ----------
    selection : InputBoolean
        Selection of atoms to apply this node to
    expand : InputBoolean
        Whether to expand the selection to the whole residue if at least one atom is selected

    Inputs
    ------
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.expand : BooleanSocket
        Whether to expand the selection to the whole residue if at least one atom is selected

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    """

    _name = "Select Res Whole"
    _asset_name = "Select Res Whole"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_res_whole"}

    class _Inputs(SocketAccessor):
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        expand: BooleanSocket
        """Whether to expand the selection to the whole residue if at least one atom is selected"""

    class _Outputs(SocketAccessor):
        selection: BooleanSocket
        """The calculated selection"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        selection: InputBoolean = False,
        expand: InputBoolean = True,
    ):
        super().__init__(**{"Selection": selection, "Expand": expand})

    def _build_group(self, tree):
        selection = tree.inputs.boolean(
            "Selection",
            False,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        expand = tree.inputs.boolean(
            "Expand",
            True,
            description="Whether to expand the selection to the whole residue if at least one atom is selected",
        )
        selection_1 = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )

        (
            expand.switch.boolean(
                selection, BooleanAny(boolean=selection, group_id=UniqueResidueID())
            )
            >> selection_1
        )


ASSET = SelectResWhole

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
