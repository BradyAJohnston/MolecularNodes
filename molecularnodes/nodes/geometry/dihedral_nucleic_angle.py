# Node-group asset "Dihedral Nucleic Angle" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .atom_name import AtomName
from .dihedral_angle import DihedralAngle
from .find_bonded_atom import FindBondedAtom
from .residue_mask import ResidueMask


class DihedralNucleicAngle(AssetGeometryGroup):
    """
    Dihedral Nucleic Angle

    Outputs
    -------
    o.angle : FloatSocket
        The angle between the vectors AB and CD, when made perpendicular to BC.
    o.up : VectorSocket
        The vector BA when made perpendicular to the axis BC
    o.axis : VectorSocket
        The axis vector BC
    """

    _name = "Dihedral Nucleic Angle"
    _asset_name = "Dihedral Nucleic Angle"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        angle: FloatSocket
        """The angle between the vectors AB and CD, when made perpendicular to BC."""
        up: VectorSocket
        """The vector BA when made perpendicular to the axis BC"""
        axis: VectorSocket
        """The axis vector BC"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        angle = tree.outputs.float(
            "Angle",
            description="The angle between the vectors AB and CD, when made perpendicular to BC.",
            subtype="ANGLE",
        )
        up = tree.outputs.vector(
            "Up", description="The vector BA when made perpendicular to  the axis BC"
        )
        axis = tree.outputs.vector("Axis", description="The axis vector BC")

        group = ResidueMask(atom_name=54)
        group_1 = ResidueMask(atom_name=50)
        group_2 = ResidueMask(atom_name=55)
        group_3 = ResidueMask(atom_name=53)
        group_4 = ResidueMask(atom_name=57)
        integer_math = AtomName().o.atom_name - 50
        index_switch = g.IndexSwitch.vector(
            integer_math,
            (
                FindBondedAtom(atom_name="O5'").o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group.o.position,
                group_2.o.position,
                group_4.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                FindBondedAtom(atom_name="P", distance=1).o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
            ),
        )
        index_switch_1 = g.IndexSwitch.vector(
            integer_math,
            (
                FindBondedAtom(atom_name="O3'", distance=1).o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_1.o.position,
                group_3.o.position,
                group.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_4.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
            ),
        )
        index_switch_2 = g.IndexSwitch.vector(
            integer_math,
            (
                FindBondedAtom(atom_name="C3'", distance=1).o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                FindBondedAtom(atom_name="O3'").o.position,
                group_1.o.position,
                group_3.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                group_2.o.position,
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
            ),
        )
        group_5 = DihedralAngle(
            a=index_switch, b=g.Position(), c=index_switch_1, d=index_switch_2
        )
        group_5.o.angle * -1.0 >> angle

        group_5.o.ba_bc >> up
        group_5.o.bc >> axis


ASSET = DihedralNucleicAngle

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
