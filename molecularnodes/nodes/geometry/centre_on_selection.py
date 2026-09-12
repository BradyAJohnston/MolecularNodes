# Node-group asset "Centre on Selection" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputGeometry, InputInteger
from .centroid import Centroid


class CentreOnSelection(AssetGeometryGroup):
    """
    Centre on Selection

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection within the groups to calculate the centroid for, which then affects all other points in the group
    group_id : InputInteger
        Optionally centre the points separately for each group based on the `Group ID` input

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection within the groups to calculate the centroid for, which then affects all other points in the group
    i.group_id : IntegerSocket
        Optionally centre the points separately for each group based on the `Group ID` input

    Outputs
    -------
    o.atoms : GeometrySocket
        Atoms that have been moved to object origin based on their group's calculated centroid
    o.offset : VectorSocket
        The calculated vector offset that was applied to the points
    """

    _name = "Centre on Selection"
    _asset_name = "Centre on Selection"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.centre_on_selection"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection within the groups to calculate the centroid for, which then affects all other points in the group"""
        group_id: IntegerSocket
        """Optionally centre the points separately for each group based on the `Group ID` input"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """Atoms that have been moved to object origin based on their group's calculated centroid"""
        offset: VectorSocket
        """The calculated vector offset that was applied to the points"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        group_id: InputInteger = 0,
    ):
        super().__init__(
            **{"Atoms": atoms, "Selection": selection, "Group ID": group_id}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms",
            description="Atomic geometry that contains vertices and edges",
            hide_value=True,
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection within the groups to calculate the centroid for, which then affects all other points in the group",
            hide_value=True,
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="Optionally centre the points separately for each group based on the `Group ID` input",
            hide_value=True,
        )
        atoms_1 = tree.outputs.geometry(
            "Atoms",
            description="Atoms that have been moved to object origin based on their group's calculated centroid",
        )
        offset = tree.outputs.vector(
            "Offset",
            description="The calculated vector offset that was applied to the points",
        )

        capture = g.CaptureAttribute.point(geometry=atoms)
        vector = capture.items.vector(
            "Vector", Centroid(selection=selection, group_id=group_id)
        )
        vector_math = vector.output * -1.0
        capture.o.geometry >> g.SetPosition(offset=vector_math) >> atoms_1

        vector_math >> offset


ASSET = CentreOnSelection

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
