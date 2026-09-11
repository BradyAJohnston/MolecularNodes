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
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
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

        group = MN_chi_atom_names()
        group_1 = AtomName()
        group_2 = ResidueMask(atom_name=group.o.x3)
        group_3 = ResidueMask(atom_name=group.o.x4)
        group_4 = MenuResidueMask(atom_name="CA")
        group_5 = MenuResidueMask(atom_name="CB")
        group_6 = MenuResidueMask(atom_name="CG")
        group_7 = MenuResidueMask(atom_name="CD")
        index_switch = g.IndexSwitch.vector(
            group_1,
            (
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                ResidueMask(atom_name=group.o.x1).o.position,
                ResidueMask(atom_name=group.o.x2).o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_2.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_2.o.position,
                group_3.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_3.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                ResidueMask(atom_name=group.o.x5).o.position,
            ),
        )
        index_switch_1 = g.IndexSwitch.vector(
            group_1,
            (
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_4.o.position,
                group_5.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_6.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_6.o.position,
                group_7.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_7.o.position,
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
        index_switch_2 = g.IndexSwitch.vector(
            group_1,
            (
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                MenuResidueMask().o.position,
                group_4.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_5.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_5.o.position,
                group_6.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_6.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_7.o.position,
            ),
        )
        group_8 = DihedralAngle(
            a=index_switch, b=g.Position(), c=index_switch_1, d=index_switch_2
        )
        (g.EdgesOfVertex().o.total > 1).switch.float(
            true=group_8.o.angle
        ) * -1.0 >> angle

        group_8.o.ba_bc >> up
        group_8.o.bc >> axis


ASSET = DihedralChiAngle

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
