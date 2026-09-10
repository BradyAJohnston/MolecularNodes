# Node-group asset "Residue Dihedral Angle" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputInteger
from .dihedral_angle import DihedralAngle
from .residue_mask import ResidueMask


class ResidueDihedralAngle(AssetGeometryGroup):
    """
    Residue Dihedral Angle

    Parameters
    ----------
    a : InputInteger
        Atom to pick from the group
    c : InputInteger
        Atom to pick from the group
    d : InputInteger
        Atom to pick from the group

    Inputs
    ------
    i.a : IntegerSocket
        Atom to pick from the group
    i.c : IntegerSocket
        Atom to pick from the group
    i.d : IntegerSocket
        Atom to pick from the group

    Outputs
    -------
    o.value : FloatSocket
        Value
    o.ba_bc : VectorSocket
        The vector BA when made perpendicular to the axis BC
    o.cd_bc : VectorSocket
        The Vector CD when makde perpendicular to the axis BC
    o.bc : VectorSocket
        The axis vector BC
    """

    _name = "Residue Dihedral Angle"
    _asset_name = "Residue Dihedral Angle"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")

    class _Inputs(SocketAccessor):
        a: IntegerSocket
        """Atom to pick from the group"""
        c: IntegerSocket
        """Atom to pick from the group"""
        d: IntegerSocket
        """Atom to pick from the group"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """Value"""
        ba_bc: VectorSocket
        """The vector BA when made perpendicular to the axis BC"""
        cd_bc: VectorSocket
        """The Vector CD when makde perpendicular to the axis BC"""
        bc: VectorSocket
        """The axis vector BC"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        a: InputInteger = 6,
        c: InputInteger = 2,
        d: InputInteger = 3,
    ):
        super().__init__(**{"A": a, "C": c, "D": d})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        a = tree.inputs.integer(
            "A", 6, description="Atom to pick from the group", min_value=2
        )
        c_ = tree.inputs.integer(
            "C", 2, description="Atom to pick from the group", min_value=2
        )
        d = tree.inputs.integer(
            "D", 3, description="Atom to pick from the group", min_value=2
        )
        value = tree.outputs.float("Value")
        ba_bc = tree.outputs.vector(
            "BA⟂(BC)",
            description="The vector BA when made perpendicular to  the axis BC",
        )
        cd_bc = tree.outputs.vector(
            "CD⟂(BC)",
            description="The Vector CD when makde perpendicular to the axis BC",
        )
        bc = tree.outputs.vector("BC", description="The axis vector BC")

        group = DihedralAngle(
            a=ResidueMask(atom_name=a).o.position,
            b=g.Position(),
            c=ResidueMask(atom_name=c_).o.position,
            d=ResidueMask(atom_name=d).o.position,
        )

        group >> value
        group.o.ba_bc >> ba_bc
        group.o.cd_bc >> cd_bc
        group.o.bc >> bc


ASSET = ResidueDihedralAngle

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
