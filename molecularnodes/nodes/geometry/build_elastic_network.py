# Node-group asset "Build Elastic Network" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry, InputMenu
from ._shared.sample_atomic_attributes import SampleAtomicAttributes
from .is_alpha_carbon import IsAlphaCarbon
from .plexus import Plexus


class BuildElasticNetwork(AssetGeometryGroup):
    """
    Build Elastic Network

    Parameters
    ----------
    atoms : InputGeometry
        Atomic mesh which inculdes vertices and edges to generate an elastic network on
    selection : InputBoolean
        Atoms which aren't selected will not be included in the final output mesh
    menu : InputMenu | Literal["Alpha Carbon", "All Atom"]
        Output a network between just alpha carbons, or one that includes all atoms
    alpha_carbon : InputFloat
        Scale the distance for edge creation between just alpha carbon atoms, ignoring all other atoms.
    all_atom : InputFloat
        Scale the distance for edge creation, looking for all selected atoms

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic mesh which inculdes vertices and edges to generate an elastic network on
    i.selection : BooleanSocket
        Atoms which aren't selected will not be included in the final output mesh
    i.menu : MenuSocket
        Output a network between just alpha carbons, or one that includes all atoms
    i.alpha_carbon : FloatSocket
        Scale the distance for edge creation between just alpha carbon atoms, ignoring all other atoms.
    i.all_atom : FloatSocket
        Scale the distance for edge creation, looking for all selected atoms

    Outputs
    -------
    o.mesh : GeometrySocket
        The generated elastic network. Edges are formed between atoms within the cutoff distance for use as constraints in a simulation
    """

    _name = "Build Elastic Network"
    _asset_name = "Build Elastic Network"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic mesh which inculdes vertices and edges to generate an elastic network on"""
        selection: BooleanSocket
        """Atoms which aren't selected will not be included in the final output mesh"""
        menu: MenuSocket
        """Output a network between just alpha carbons, or one that includes all atoms"""
        alpha_carbon: FloatSocket
        """Scale the distance for edge creation between just alpha carbon atoms, ignoring all other atoms."""
        all_atom: FloatSocket
        """Scale the distance for edge creation, looking for all selected atoms"""

    class _Outputs(SocketAccessor):
        mesh: GeometrySocket
        """The generated elastic network. Edges are formed between atoms within the cutoff distance for use as constraints in a simulation"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        menu: InputMenu | Literal["Alpha Carbon", "All Atom"] = "Alpha Carbon",
        alpha_carbon: InputFloat = 4.0,
        all_atom: InputFloat = 2.0,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Menu": menu,
                "Alpha Carbon": alpha_carbon,
                "All Atom": all_atom,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms",
            description="Atomic mesh which inculdes vertices and edges to generate an elastic network on",
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Atoms which aren't selected will not be included in the final output mesh",
            hide_value=True,
        )
        menu = tree.inputs.menu(
            "Menu",
            description="Output a network between just alpha carbons, or one that includes all atoms",
            optional_label=True,
        )
        alpha_carbon = tree.inputs.float(
            "Alpha Carbon",
            4.0,
            description="Scale the distance for edge creation between just alpha carbon atoms, ignoring all other atoms.",
            min_value=0.0,
            max_value=10_000.0,
        )
        all_atom = tree.inputs.float(
            "All Atom",
            2.0,
            description="Scale the distance for edge creation, looking for all selected atoms",
            min_value=0.0,
            max_value=10_000.0,
        )
        mesh = tree.outputs.geometry(
            "Mesh",
            description="The generated elastic network. Edges are formed between atoms within the cutoff distance for use as constraints in a simulation",
        )

        separate_geometry = g.SeparateGeometry.point(atoms, selection)
        group = Plexus(
            points=g.SeparateGeometry.point(
                separate_geometry.o.selection, IsAlphaCarbon().o.selection
            ).o.selection,
            distance=alpha_carbon,
            radius=0.0,
        )
        group_1 = Plexus(
            points=g.SeparateGeometry.point(separate_geometry.o.selection).o.selection,
            distance=all_atom,
            radius=0.0,
        )
        merge_by_distance = g.MenuSwitch.geometry(
            menu,
            {
                "Alpha Carbon": group,
                "All Atom": g.JoinGeometry(geometry=(group, group_1)),
            },
        ) >> g.MergeByDistance(distance=0.0001)
        capture = g.CaptureAttribute.point(geometry=merge_by_distance)
        index = capture.items.integer(
            "Index", g.SampleNearest.point(separate_geometry.o.selection)
        )
        (
            SampleAtomicAttributes(
                atoms=capture.o.geometry,
                sample_atoms=separate_geometry.o.selection,
                index=index.output,
            )
            >> g.SortElements.point(sort_weight=index.output)
            >> mesh
        )

        menu.default_value = "Alpha Carbon"


ASSET = BuildElasticNetwork

ASSET_METADATA = {
    "catalog_id": "c2c958af-5095-4fc2-884d-709bba965fc4",
}
