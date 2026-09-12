# Node-group asset "Find Bonds" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from ._shared.sample_atomic_attributes import SampleAtomicAttributes
from .evaluate_on_atoms import EvaluateOnAtoms
from .plexus import Plexus
from .vdw_radii import VDWRadii


class FindBonds(AssetGeometryGroup):
    """
    Find Bonds

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    scale : InputFloat
        Scale the VDW radii of the atoms when searching for bonds

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.scale : FloatSocket
        Scale the VDW radii of the atoms when searching for bonds

    Outputs
    -------
    o.atoms : GeometrySocket
        Atoms
    """

    _name = "Find Bonds"
    _asset_name = "Find Bonds"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.topology_find_bonds"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        scale: FloatSocket
        """Scale the VDW radii of the atoms when searching for bonds"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """Atoms"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        scale: InputFloat = 1.0,
    ):
        super().__init__(**{"Atoms": atoms, "Selection": selection, "Scale": scale})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        scale = tree.inputs.float(
            "Scale",
            1.0,
            description="Scale the VDW radii of the atoms when searching for bonds",
            min_value=0.0,
            max_value=10_000.0,
        )
        atoms_1 = tree.outputs.geometry("Atoms")

        closure_zone = g.ClosureZone()
        atoms_2 = closure_zone.inputs.geometry("Atoms")
        geometry = closure_zone.outputs.geometry("Geometry")
        separate_geometry = atoms_2 >> g.SeparateGeometry.point(selection=selection)
        group = Plexus(
            points=separate_geometry.o.selection,
            distance=scale,
            radius=VDWRadii().o.vdw_radii * 0.58,
        )
        (
            SampleAtomicAttributes(
                atoms=group, sample_atoms=separate_geometry.o.selection
            )
            >> geometry
        )
        EvaluateOnAtoms(geometry=atoms, closure=closure_zone.closure) >> atoms_1


ASSET = FindBonds

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
