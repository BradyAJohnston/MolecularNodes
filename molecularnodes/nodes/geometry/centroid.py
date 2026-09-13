# Node-group asset "Centroid" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputInteger, InputVector


class Centroid(AssetGeometryGroup):
    """
    Centroid

    Parameters
    ----------
    position : InputVector
        The `Position` vector to use for the centroid calculation
    selection : InputBoolean
        Selected points contribute to the computation of the centroid, unselected points do not contribute but still return the centroid for their `Group ID`
    group_id : InputInteger
        Compute the centroid on for each unique `Group ID`

    Inputs
    ------
    i.position : VectorSocket
        The `Position` vector to use for the centroid calculation
    i.selection : BooleanSocket
        Selected points contribute to the computation of the centroid, unselected points do not contribute but still return the centroid for their `Group ID`
    i.group_id : IntegerSocket
        Compute the centroid on for each unique `Group ID`

    Outputs
    -------
    o.centroid : VectorSocket
        The computed average vector for each `Group ID`
    """

    _name = "Centroid"
    _asset_name = "Centroid"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.centroid"}

    class _Inputs(SocketAccessor):
        position: VectorSocket
        """The `Position` vector to use for the centroid calculation"""
        selection: BooleanSocket
        """Selected points contribute to the computation of the centroid, unselected points do not contribute but still return the centroid for their `Group ID`"""
        group_id: IntegerSocket
        """Compute the centroid on for each unique `Group ID`"""

    class _Outputs(SocketAccessor):
        centroid: VectorSocket
        """The computed average vector for each `Group ID`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        position: InputVector = None,
        selection: InputBoolean = True,
        group_id: InputInteger = 0,
    ):
        super().__init__(
            **{"Position": position, "Selection": selection, "Group ID": group_id}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        position = tree.inputs.vector(
            "Position",
            (0.0, 0.0, 0.0),
            description="The `Position` vector to use for the centroid calculation",
            hide_value=True,
            default_input="POSITION",
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selected points contribute to the computation of the centroid, unselected points do not contribute but still return the centroid for their `Group ID`",
            hide_value=True,
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="Compute the centroid on for each unique `Group ID`",
            hide_value=True,
        )
        centroid = tree.outputs.vector(
            "Centroid", description="The computed average vector for each `Group ID`"
        )

        (
            selection.switch.vector((0.0, 0.0, 0.0), position).point.total(group_id)
            / g.AccumulateField.point.integer(selection, group_id).o.total
            >> centroid
        )


ASSET = Centroid

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
