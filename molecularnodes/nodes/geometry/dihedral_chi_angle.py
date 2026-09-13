# Node-group asset "Dihedral Chi Angle" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from ._shared.mn_chi_atom_names import MN_chi_atom_names
from .atom_name import AtomName
from .dihedral_angle import DihedralAngle
from .menu_residue_mask import MenuResidueMask
from .residue_mask import ResidueMask


class DihedralChiAngle(AssetGeometryGroup):
    """
    Dihedral Chi Angle

    Outputs
    -------
    o.angle : FloatSocket
        Angle
    o.up : VectorSocket
        Up
    o.axis : VectorSocket
        Axis
    """

    _name = "Dihedral Chi Angle"
    _asset_name = "Dihedral Chi Angle"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        angle: FloatSocket
        """Angle"""
        up: VectorSocket
        """Up"""
        axis: VectorSocket
        """Axis"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        angle = tree.outputs.float("Angle", subtype="ANGLE")
        up = tree.outputs.vector("Up")
        axis = tree.outputs.vector("Axis")

        mn_chi_atom_names = MN_chi_atom_names()
        atom_name = AtomName()
        menu_residue_mask = MenuResidueMask(atom_name="CA")
        menu_residue_mask_1 = MenuResidueMask(atom_name="CB")
        menu_residue_mask_2 = MenuResidueMask(atom_name="CG")
        menu_residue_mask_3 = MenuResidueMask(atom_name="CD")
        index_switch = g.IndexSwitch.vector(
            atom_name,
            (
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                menu_residue_mask.o.position,
                menu_residue_mask_1.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                menu_residue_mask_2.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                menu_residue_mask_2.o.position,
                menu_residue_mask_3.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                menu_residue_mask_3.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                MenuResidueMask(atom_name="NE").o.position,
            ),
        )
        index_switch_1 = g.IndexSwitch.vector(
            atom_name,
            (
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                MenuResidueMask().o.position,
                menu_residue_mask.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                menu_residue_mask_1.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                menu_residue_mask_1.o.position,
                menu_residue_mask_2.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                menu_residue_mask_2.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                menu_residue_mask_3.o.position,
            ),
        )
        residue_mask = ResidueMask(atom_name=mn_chi_atom_names.o.x3)
        residue_mask_1 = ResidueMask(atom_name=mn_chi_atom_names.o.x4)
        index_switch_2 = g.IndexSwitch.vector(
            atom_name,
            (
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                ResidueMask(atom_name=mn_chi_atom_names.o.x1).o.position,
                ResidueMask(atom_name=mn_chi_atom_names.o.x2).o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                residue_mask.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                residue_mask.o.position,
                residue_mask_1.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                residue_mask_1.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                ResidueMask(atom_name=mn_chi_atom_names.o.x5).o.position,
            ),
        )
        dihedral_angle = DihedralAngle(
            a=index_switch_2, b=g.Position(), c=index_switch, d=index_switch_1
        )
        (
            (g.EdgesOfVertex().o.total > 1).switch.float(true=dihedral_angle.o.angle)
            * -1.0
            >> angle
        )

        dihedral_angle.o.ba_bc >> up
        dihedral_angle.o.bc >> axis


ASSET = DihedralChiAngle

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
