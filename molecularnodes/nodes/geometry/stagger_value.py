# Node-group asset "Stagger Value" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
)
from nodebpy.types import Default, InputFloat, InputInteger


class StaggerValue(AssetGeometryGroup):
    """
    Stagger Value

    Parameters
    ----------
    value : InputFloat
        Animated value that interpolates from min to max over frames
    width : InputFloat
        How many `ID`s are staggered at the same time during the animation.
    id : InputInteger
        The `ID` over which we stagger the animation. When unconnected: The "id" attribute if available, otherwise the index.
    group_id : InputInteger
        Compute the staggered animation individually for each `Group ID`

    Inputs
    ------
    i.value : FloatSocket
        Animated value that interpolates from min to max over frames
    i.width : FloatSocket
        How many `ID`s are staggered at the same time during the animation.
    i.id : IntegerSocket
        The `ID` over which we stagger the animation.
    i.group_id : IntegerSocket
        Compute the staggered animation individually for each `Group ID`

    Outputs
    -------
    o.value : FloatSocket
        Value
    """

    _name = "Stagger Value"
    _asset_name = "Stagger Value"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        value: FloatSocket
        """Animated value that interpolates from min to max over frames"""
        width: FloatSocket
        """How many `ID`s are staggered at the same time during the animation."""
        id: IntegerSocket
        """The `ID` over which we stagger the animation."""
        group_id: IntegerSocket
        """Compute the staggered animation individually for each `Group ID`"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """Value"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputFloat = 0.0,
        width: InputFloat = 3.0,
        id: InputInteger = Default.ID_OR_INDEX,
        group_id: InputInteger = 0,
    ):
        super().__init__(
            **{"Value": value, "Width": width, "ID": id, "Group ID": group_id}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        value = tree.inputs.float(
            "Value",
            0.0,
            description="Animated value that interpolates from min to max over frames",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        width = tree.inputs.float(
            "Width",
            3.0,
            description="How many `ID`s are staggered at the same time during the animation.",
            min_value=0.0,
            max_value=10_000.0,
        )
        id = tree.inputs.integer(
            "ID",
            0,
            description="The `ID` over which we stagger the animation.",
            default_input="ID_OR_INDEX",
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="Compute the staggered animation individually for each `Group ID`",
            hide_value=True,
        )
        value_1 = tree.outputs.float("Value")

        field_min_max = g.FieldMinAndMax.point.integer(id, group_id)
        math_1 = g.Math.subtract(id, field_min_max.o.min)
        math_2 = (
            g.Math.subtract(field_min_max.o.max, field_min_max.o.min).o.value + width
        )
        math_3 = math_1.o.value / math_2
        (
            (width > 0.0).switch.float(
                value >= math_3,
                value.map_range(math_3, (math_1.o.value + width) / math_2),
            )
            >> value_1
        )


ASSET = StaggerValue

ASSET_METADATA = {
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
