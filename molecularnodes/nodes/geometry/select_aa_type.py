# Node-group asset 'Select AA Type' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean
from .boolean_andor import BooleanAndOr
from .residue_name import ResidueName


class SelectAaType(AssetGeometryGroup):
    """
    Select AA Type

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection

    Outputs
    -------
    o.is_polar : BooleanSocket
        is_polar
    o.is_apolar : BooleanSocket
        is_apolar
    o.is_acidic : BooleanSocket
        is_acidic
    o.is_basic : BooleanSocket
        is_basic
    o.is_aromatic : BooleanSocket
        is_aromatic
    """

    _name = "Select AA Type"
    _asset_name = "Select AA Type"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_aa_type"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""

    class _Outputs(SocketAccessor):
        is_polar: BooleanSocket
        is_apolar: BooleanSocket
        is_acidic: BooleanSocket
        is_basic: BooleanSocket
        is_aromatic: BooleanSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        and_: InputBoolean = True,
        or_: InputBoolean = False,
    ):
        super().__init__(**{"And": and_, "Or": or_})

    def _build_group(self, tree):
        and_ = tree.inputs.boolean(
            "And",
            True,
            description="The resulting selection must overlap with this input selection",
            hide_value=True,
        )
        or_ = tree.inputs.boolean(
            "Or",
            False,
            description="The resulting selection can be calculated from this node or be from this input selection",
            hide_value=True,
        )
        is_polar = tree.outputs.boolean("is_polar")
        is_apolar = tree.outputs.boolean("is_apolar")
        is_acidic = tree.outputs.boolean("is_acidic")
        is_basic = tree.outputs.boolean("is_basic")
        is_aromatic = tree.outputs.boolean("is_aromatic")

        boolean = g.Boolean(boolean=True)
        group = ResidueName()

        # Polar: ASN(2), CYS(4), GLN(6), HIS(8), SER(15), THR(16)
        index_switch_polar = g.IndexSwitch.boolean(
            group,
            (
                False, False, boolean, False, boolean, False,
                boolean, False, boolean, False, False, False,
                False, False, False, boolean, boolean, False,
                False, False,
            ) + (False,) * 24,  # pad to 44 for DNA/RNA
        )
        BooleanAndOr(and_=and_, or_=or_, boolean=index_switch_polar) >> is_polar

        # Apolar: ALA(0), GLY(7), ILE(9), LEU(10), MET(12), PRO(14), VAL(19)
        index_switch_apolar = g.IndexSwitch.boolean(
            group,
            (
                boolean, False, False, False, False, False,
                False, boolean, False, boolean, boolean, False,
                boolean, False, boolean, False, False, False,
                False, boolean,
            ) + (False,) * 24,
        )
        BooleanAndOr(and_=and_, or_=or_, boolean=index_switch_apolar) >> is_apolar

        # Acidic: ASP(3), GLU(5)
        index_switch_acidic = g.IndexSwitch.boolean(
            group,
            (
                False, False, False, boolean, False, boolean,
                False, False, False, False, False, False,
                False, False, False, False, False, False,
                False, False,
            ) + (False,) * 24,
        )
        BooleanAndOr(and_=and_, or_=or_, boolean=index_switch_acidic) >> is_acidic

        # Basic: ARG(1), LYS(11)
        index_switch_basic = g.IndexSwitch.boolean(
            group,
            (
                False, boolean, False, False, False, False,
                False, False, False, False, False, boolean,
                False, False, False, False, False, False,
                False, False,
            ) + (False,) * 24,
        )
        BooleanAndOr(and_=and_, or_=or_, boolean=index_switch_basic) >> is_basic

        # Aromatic: PHE(13), TRP(17), TYR(18)
        index_switch_aromatic = g.IndexSwitch.boolean(
            group,
            (
                False, False, False, False, False, False,
                False, False, False, False, False, False,
                False, boolean, False, False, False, boolean,
                boolean, False,
            ) + (False,) * 24,
        )
        BooleanAndOr(and_=and_, or_=or_, boolean=index_switch_aromatic) >> is_aromatic


ASSET = SelectAaType

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
