# Node-group asset "Topology DSSP" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    CustomGeometryGroup,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputGeometry
from ._shared.mn_topo_assign_backbone import MN_topo_assign_backbone
from .backbone_c import BackboneC
from .backbone_n import BackboneN
from .backbone_nh import BackboneNH
from .backbone_o import BackboneO
from .boolean_run_fill import BooleanRunFill
from .boolean_run_trim import BooleanRunTrim
from .evaluate_on_atoms import EvaluateOnAtoms
from .integer_distance import IntegerDistance
from .is_alpha_carbon import IsAlphaCarbon
from .menu_residue_mask import MenuResidueMask
from .offset_boolean import OffsetBoolean
from .offset_index import OffsetIndex
from .offset_vector import OffsetVector
from .secondary_structure import SecondaryStructure
from .ures_id import UResID
from .visualize_relative_atoms import VisualizeRelativeAtoms
from .world_to_angstrom import WorldToAngstrom


class RecipAngDis(CustomGeometryGroup):
    _name = "Recip. Ang. Dis"
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.recip_ang_dis"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        a = tree.inputs.vector(
            "A", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0
        )
        b = tree.inputs.vector(
            "B", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0
        )
        value = tree.outputs.float("Value")

        1.0 / WorldToAngstrom(world=a.distance(b)) >> value


class HBondEnergy(CustomGeometryGroup):
    _name = "HBond Energy"
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.hbond_energy"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        o = tree.inputs.vector("O", (0.0, 0.0, 0.0))
        c_ = tree.inputs.vector("C", (0.0, 0.0, 0.0))
        n = tree.inputs.vector("N", (0.0, 0.0, 0.0))
        h = tree.inputs.vector("H", (0.0, 0.0, 0.0))
        is_bonded = tree.outputs.boolean("Is Bonded")
        bond_energy = tree.outputs.float("Bond Energy")
        bond_vector = tree.outputs.vector("Bond Vector")

        with g.Frame("1 / r(ON)"):
            group = RecipAngDis(A=o, B=n)
        with g.Frame("1 / r(CH)"):
            group_1 = RecipAngDis(A=c_, B=h)
        with g.Frame("1 / r(OH)"):
            group_2 = RecipAngDis(A=o, B=h)
        with g.Frame("1 / r(CN)"):
            group_3 = RecipAngDis(A=c_, B=n)
        math_1 = (group.o.value + group_1 - group_2 - group_3) * -1.0
        math_1.node.mute = True
        math_2 = math_1 * 0.084 * 332.0
        (math_2 < -0.5) >> is_bonded
        o - h >> bond_vector

        math_2 >> bond_energy


class CheckHBond(CustomGeometryGroup):
    _name = ".Check HBond"
    _color_tag = "CONVERTER"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        co_index = tree.inputs.integer(
            "CO Index", 0, min_value=0, default_input="INDEX"
        )
        nh_index = tree.inputs.integer(
            "NH Index", 0, min_value=0, default_input="INDEX"
        )
        distance = tree.inputs.integer("Distance", 2)
        is_bonded = tree.outputs.boolean("Is Bonded")
        bond_energy = tree.outputs.float("Bond Energy")
        bond_vector = tree.outputs.vector("Bond Vector")

        group = HBondEnergy(
            O=BackboneO(method="Read").o.o.point.at(co_index),
            C=BackboneC(method="Read").o.c.point.at(co_index),
            N=BackboneN(method="Read").o.n.point.at(nh_index),
            H=BackboneNH(menu="Read").o.nh.point.at(nh_index),
        )
        group_1 = IntegerDistance(
            a=UResID(index=co_index).o.ures_id,
            b=UResID(index=nh_index).o.ures_id,
            distance=distance,
        )
        (group.o.is_bonded & group_1.o.cutoff) >> is_bonded

        group.o.bond_energy >> bond_energy
        group.o.bond_vector >> bond_vector


class HBondBackboneCheck(CustomGeometryGroup):
    _name = ".HBond Backbone Check"
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry._hbond_backbone_check"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        with tree.inputs.panel("CO (i)"):
            co_index = tree.inputs.integer(
                "CO Index", 0, min_value=0, default_input="INDEX"
            )
            co_offset = tree.inputs.integer("CO Offset", 0)
        with tree.inputs.panel("NH (j)"):
            nh_index = tree.inputs.integer(
                "NH Index", 0, min_value=0, default_input="INDEX"
            )
            nh_offset = tree.inputs.integer("NH Offset", 0)
        is_bonded = tree.outputs.boolean("Is Bonded")
        bond_energy = tree.outputs.float("Bond Energy")
        h_o = tree.outputs.vector("H->O")

        group = CheckHBond(
            **{"CO Index": co_index + co_offset, "NH Index": nh_index + nh_offset},
            Distance=1,
        )

        group >> is_bonded
        group.o.bond_energy >> bond_energy
        group.o.bond_vector >> h_o


class MN_topo_calc_helix(CustomGeometryGroup):
    _name = ".MN_topo_calc_helix"
    _tree_properties = {"node_tool_idname": "geometry._mn_topo_calc_helix"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        ca_mesh = tree.inputs.geometry("CA Mesh")
        is_helix = tree.outputs.boolean("Is Helix")
        instances = tree.outputs.geometry("Instances")

        group = BackboneNH(menu="Read")
        capture = g.CaptureAttribute.point(geometry=ca_mesh)
        n_3_helix = capture.items.boolean(
            "3-helix", HBondBackboneCheck(**{"NH Offset": 3}).o.is_bonded
        )
        n_4_helix = capture.items.boolean(
            "4-helix", HBondBackboneCheck(**{"NH Offset": 4}).o.is_bonded
        )
        n_5_helix = capture.items.boolean(
            "5-helix", HBondBackboneCheck(**{"NH Offset": 6}).o.is_bonded
        )
        with g.Frame():
            boolean_math = n_4_helix.output & OffsetBoolean(
                boolean=n_4_helix.output, offset=1
            )
            boolean_math_1 = (
                boolean_math
                | OffsetBoolean(boolean=boolean_math, offset=-1)
                | OffsetBoolean(boolean=boolean_math, offset=-2)
            )
            boolean_math_2 = (
                boolean_math_1
                | OffsetBoolean(boolean=boolean_math, offset=-3)
                | OffsetBoolean(boolean=boolean_math, offset=-4)
            )
            boolean_math_3 = boolean_math_2 | OffsetBoolean(
                boolean=boolean_math, offset=-5
            )
        with g.Frame():
            boolean_math_4 = n_5_helix.output & OffsetBoolean(
                boolean=n_5_helix.output, offset=1
            )
            boolean_math_5 = (
                boolean_math_4
                | OffsetBoolean(boolean=boolean_math_4, offset=-1)
                | OffsetBoolean(boolean=boolean_math_4, offset=-2)
            )
            boolean_math_6 = (
                boolean_math_5
                | OffsetBoolean(boolean=boolean_math_4, offset=-3)
                | OffsetBoolean(boolean=boolean_math_4, offset=-4)
            )
            boolean_math_7 = (
                boolean_math_6
                | OffsetBoolean(boolean=boolean_math_4, offset=-5)
                | OffsetBoolean(boolean=boolean_math_4, offset=-6)
            )
        with g.Frame():
            boolean_math_8 = n_3_helix.output & OffsetBoolean(
                boolean=n_3_helix.output, offset=1
            )
            boolean_math_9 = (
                boolean_math_8
                | OffsetBoolean(boolean=boolean_math_8, offset=-1)
                | OffsetBoolean(boolean=boolean_math_8, offset=-2)
            )
            boolean_math_10 = (
                boolean_math_9
                | OffsetBoolean(boolean=boolean_math_8, offset=-3)
                | OffsetBoolean(boolean=boolean_math_8, offset=-4)
            )
        group_1 = OffsetIndex(offset=4)
        _group_2 = OffsetVector(vector=group, index=group_1)
        sample_index = g.SampleIndex(
            geometry=capture.o.geometry,
            value=boolean_math_10 | boolean_math_3 | boolean_math_7,
            index=g.Index(),
            data_type="BOOLEAN",
        )
        (
            VisualizeRelativeAtoms(
                atoms=capture.o.geometry,
                selection=n_4_helix.output,
                position=BackboneO(method="Read"),
                target_index=group_1,
                target_position=group,
            )
            >> instances
        )

        sample_index >> is_helix


class MN_topo_calc_sheet(CustomGeometryGroup):
    _name = ".MN_topo_calc_sheet"
    _tree_properties = {"node_tool_idname": "geometry._mn_topo_calc_sheet"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        ca_mesh = tree.inputs.geometry("CA Mesh")
        is_sheet = tree.outputs.boolean("Is Sheet")
        instances = tree.outputs.geometry("Instances")

        with g.Frame("Find Residues that we might be HBonding to"):
            group = BackboneNH(menu="Read")
            group_1 = BackboneO(method="Read")
            capture = g.CaptureAttribute.point(geometry=ca_mesh)
            o_nh = capture.items.integer(
                "O -> NH",
                g.SampleNearest.point(
                    g.SetPosition(geometry=ca_mesh, position=group), group_1
                ),
            )
            nh_o = capture.items.integer(
                "NH -> O",
                g.SampleNearest.point(
                    g.SetPosition(geometry=ca_mesh, position=group_1), group
                ),
            )
        with g.Frame("Check if they are actually bonded to to the relevant atom"):
            capture_1 = g.CaptureAttribute.point(geometry=capture.o.geometry)
            co_nh = capture_1.items.boolean(
                "CO:NH", CheckHBond(**{"NH Index": o_nh.output}, Distance=3).o.is_bonded
            )
            nh_co = capture_1.items.boolean(
                "NH:CO", CheckHBond(**{"CO Index": nh_o.output}, Distance=4).o.is_bonded
            )
        with g.Frame("Debug arrows for HBonds"):
            value = g.Value(1.0)
            group_2 = VisualizeRelativeAtoms(
                atoms=capture_1.o.geometry,
                selection=co_nh.output,
                scale=value,
                position=BackboneO(method="Read"),
                target_index=o_nh.output,
                target_position=BackboneNH(menu="Read"),
            )
            group_3 = VisualizeRelativeAtoms(
                atoms=capture_1.o.geometry,
                selection=nh_co.output,
                scale=value,
                position=BackboneNH(menu="Read"),
                target_index=nh_o.output,
                target_position=BackboneO(method="Read"),
            )
            join_geometry = g.JoinGeometry(geometry=(group_2, group_3))
        with g.Frame("Not 100% correct but best I can do without Lists"):
            boolean_math = OffsetBoolean(
                boolean=co_nh.output, offset=-1
            ).o.boolean & OffsetBoolean(boolean=nh_co.output, offset=1)
            group_4 = BooleanRunFill(
                boolean=co_nh.output
                & OffsetBoolean(boolean=nh_co.output, index=o_nh.output)
                | boolean_math,
                fill_size=2,
            )
        capture_2 = g.CaptureAttribute.point(geometry=capture_1.o.geometry)
        boolean = capture_2.items.boolean("Boolean", group_4)
        (
            capture_2.o.geometry
            >> g.SampleIndex(value=boolean.output, index=g.Index(), data_type="BOOLEAN")
            >> is_sheet
        )

        join_geometry >> instances


class TopologyDSSP(AssetGeometryGroup):
    """
    Calculate the secondary structure attributes for the protein chains, based on the 1983 Kabsch algorithm

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
        The input `Atoms` with updated `sec_struct` attributes based on the DSSP algorithm
    """

    _name = "Topology DSSP"
    _asset_name = "Topology DSSP"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Calculate the secondary structure attributes for the protein chains, based on the 1983 Kabsch algorithm",
        "node_tool_idname": "geometry.topology_dssp",
    }

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """The input `Atoms` with updated `sec_struct` attributes based on the DSSP algorithm"""

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

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        atoms_1 = tree.outputs.geometry(
            "Atoms",
            description="The input `Atoms` with updated `sec_struct` attributes based on the DSSP algorithm",
        )

        closure_zone = g.ClosureZone()
        atoms_2 = closure_zone.inputs.geometry("Atoms")
        geometry = closure_zone.outputs.geometry("Geometry")
        group = MN_topo_assign_backbone(atoms=atoms_2)
        capture = g.CaptureAttribute.point(geometry=group.o.ca_atoms)
        is_sheet = capture.items.boolean(
            "Is Sheet", MN_topo_calc_sheet(**{"CA Mesh": group.o.ca_atoms}).o.is_sheet
        )
        is_helix = capture.items.boolean(
            "Is Helix", MN_topo_calc_helix(**{"CA Mesh": group.o.ca_atoms}).o.is_helix
        )
        switch = BooleanRunTrim(
            boolean=g.BooleanMath.subtract(is_sheet.output, is_helix.output), size=3
        ).o.boolean.switch.integer(3, 2)
        sample_index = capture.o.geometry >> g.SampleIndex(
            value=is_helix.output.switch.integer(switch, 1),
            index=group.o.sample_index,
            data_type="INT",
        )
        store_named_attribute = g.StoreNamedAttribute.point.integer(
            atoms_2, IsAlphaCarbon().o.selection, "sec_struct", sample_index
        )
        store_named_attribute_1 = g.StoreNamedAttribute.point.integer(
            store_named_attribute,
            name="sec_struct",
            value=SecondaryStructure().o.sec_struct.point.at(
                MenuResidueMask(atom_name="CA").o.index
            ),
        )
        store_named_attribute_1 >> geometry
        (
            EvaluateOnAtoms(
                geometry=atoms, closure=closure_zone.closure, result="Bundle"
            )
            >> atoms_1
        )


ASSET = TopologyDSSP

ASSET_METADATA = {
    "description": "Calculate the secondary structure attributes for the protein chains, based on the 1983 Kabsch algorithm",
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
