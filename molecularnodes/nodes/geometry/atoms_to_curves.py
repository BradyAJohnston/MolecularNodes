# Node-group asset "Atoms to Curves" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry
from .atom_id import AtomID
from .chain_id import ChainID
from .unique_chain_id import UniqueChainID


class AtomsToCurves(AssetGeometryGroup):
    """
    Atoms to Curves

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Points to maintain and use for splines, non-selected points are removed
    sort_points : InputBoolean
        Sort points first by atom_id and chain_id before splitting to curves
    cutoff : InputFloat
        Cutoff

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Points to maintain and use for splines, non-selected points are removed
    i.sort_points : BooleanSocket
        Sort points first by atom_id and chain_id before splitting to curves
    i.cutoff : FloatSocket
        Cutoff

    Outputs
    -------
    o.curves : GeometrySocket
        Curves
    """

    _name = "Atoms to Curves"
    _asset_name = "Atoms to Curves"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.atoms_to_curves"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Points to maintain and use for splines, non-selected points are removed"""
        sort_points: BooleanSocket
        """Sort points first by atom_id and chain_id before splitting to curves"""
        cutoff: FloatSocket
        """Cutoff"""

    class _Outputs(SocketAccessor):
        curves: GeometrySocket
        """Curves"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        sort_points: InputBoolean = False,
        cutoff: InputFloat = 6.0,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Sort Points": sort_points,
                "Cutoff": cutoff,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Points to maintain and use for splines, non-selected points are removed",
            hide_value=True,
        )
        sort_points = tree.inputs.boolean(
            "Sort Points",
            False,
            description="Sort points first by atom_id and chain_id before splitting to curves",
        )
        cutoff = tree.inputs.float(
            "Cutoff", 6.0, min_value=-10_000.0, max_value=10_000.0, subtype="DISTANCE"
        )
        curves = tree.outputs.geometry("Curves")

        group = ChainID()
        separate_components = g.SeparateComponents(geometry=atoms)
        join_geometry = g.JoinGeometry(
            geometry=(
                g.CurveToMesh(curve=separate_components.o.curve),
                separate_components.o.mesh,
            )
        )
        mesh_to_points = g.MeshToPoints(
            mesh=join_geometry, selection=selection, radius=0.05
        )
        sort_elements = g.SortElements.point(
            g.SortElements.point(mesh_to_points, sort_weight=group),
            group_id=group,
            sort_weight=AtomID(),
        )
        (
            sort_points.switch.geometry(mesh_to_points, sort_elements)
            >> g.PointsToCurves(curve_group_id=UniqueChainID(cutoff=cutoff))
            >> curves
        )


ASSET = AtomsToCurves

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
