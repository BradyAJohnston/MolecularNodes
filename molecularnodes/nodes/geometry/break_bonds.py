# Node-group asset "Break Bonds" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .angstrom_to_world import AngstromToWorld
from .evaluate_on_atoms import EvaluateOnAtoms


class BreakBonds(AssetGeometryGroup):
    """
    Break Bonds

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    cutoff : InputFloat
        Cutoff distance over which to remove bonds (Angstrom)

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.cutoff : FloatSocket
        Cutoff distance over which to remove bonds (Angstrom)

    Outputs
    -------
    o.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    """

    _name = "Break Bonds"
    _asset_name = "Break Bonds"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.topology_break_bonds"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        cutoff: FloatSocket
        """Cutoff distance over which to remove bonds (Angstrom)"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = False,
        cutoff: InputFloat = 2.5,
    ):
        super().__init__(**{"Atoms": atoms, "Selection": selection, "Cutoff": cutoff})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            False,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        cutoff = tree.inputs.float(
            "Cutoff",
            2.5,
            description="Cutoff distance over which to remove bonds (Angstrom)",
            min_value=0.0,
            max_value=10_000.0,
        )
        atoms_1 = tree.outputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )

        closure_zone = g.ClosureZone()
        atoms_2 = closure_zone.inputs.geometry("Atoms")
        geometry = closure_zone.outputs.geometry("Geometry")
        edge_vertices = g.EdgeVertices()
        compare = edge_vertices.o.position_1.distance(
            edge_vertices.o.position_2
        ) > AngstromToWorld(angstrom=cutoff)
        delete_geometry = atoms_2 >> g.DeleteGeometry(
            selection=selection | compare, mode="EDGE_FACE", domain="EDGE"
        )
        delete_geometry >> geometry
        EvaluateOnAtoms(geometry=atoms, closure=closure_zone.closure) >> atoms_1


ASSET = BreakBonds

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
