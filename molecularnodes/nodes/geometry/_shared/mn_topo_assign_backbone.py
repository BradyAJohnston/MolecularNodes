# Node group '.MN_topo_assign_backbone' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    GeometrySocket,
    IntegerSocket,
    SocketAccessor,
)
from nodebpy.types import InputGeometry
from ..backbone_nh import BackboneNH
from ..is_alpha_carbon import IsAlphaCarbon
from ..is_backbone import IsBackbone
from ..menu_residue_mask import MenuResidueMask


class MN_topo_assign_backbone(CustomGeometryGroup):
    """
    .MN_topo_assign_backbone

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges

    Outputs
    -------
    o.atoms : GeometrySocket
        Atoms
    o.ca_atoms : GeometrySocket
        CA Atoms
    o.sample_index : IntegerSocket
        Sample Index
    """

    _name = ".MN_topo_assign_backbone"
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry._mn_topo_assign_backbone"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """Atoms"""
        ca_atoms: GeometrySocket
        """CA Atoms"""
        sample_index: IntegerSocket
        """Sample Index"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
    ):
        super().__init__(**{"Atoms": atoms})

    def _build_group(self, tree):
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        atoms_1 = tree.outputs.geometry("Atoms")
        ca_atoms = tree.outputs.geometry("CA Atoms")
        sample_index = tree.outputs.integer("Sample Index")

        with g.Frame("Compute only on backbone atoms, but capture their idx first"):
            group = IsAlphaCarbon()
            capture = g.CaptureAttribute.point(geometry=atoms)
            selection = capture.items.boolean("Selection", group.o.selection)
            index = capture.items.integer(
                "Index", g.AccumulateField.point.integer(group.o.selection).o.trailing
            )
            separate_geometry = g.SeparateGeometry.point(
                capture.o.geometry, IsBackbone().o.selection
            )
        sample_index_1 = g.SampleIndex(
            geometry=capture.o.geometry,
            value=index.output,
            index=g.Index(),
            data_type="INT",
        )
        repeat_zone = g.RepeatZone(4)
        geometry = repeat_zone.items.geometry("Geometry", separate_geometry.o.selection)
        group_1 = MenuResidueMask(
            atom_name=g.IndexSwitch.menu(repeat_zone.iteration, ("N", "CA", "C", "O"))
        )
        join_strings = g.JoinStrings(
            (
                g.String(string="backbone"),
                g.IndexSwitch.string(repeat_zone.iteration, ("N", "CA", "C", "O")),
            ),
            delimiter="_",
        )
        store_named_attribute = g.StoreNamedAttribute.point.vector(
            geometry.current, group_1.o.is_valid, join_strings, group_1.o.position
        )
        store_named_attribute >> geometry.next
        (
            g.SeparateGeometry.point(geometry.result, selection.output)
            >> g.StoreNamedAttribute.point.vector(
                name="backbone_NH", value=BackboneNH()
            )
            >> ca_atoms
        )

        geometry.result >> atoms_1
        sample_index_1 >> sample_index
