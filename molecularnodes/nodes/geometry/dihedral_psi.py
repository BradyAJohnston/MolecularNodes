# Node-group asset "Dihedral Psi" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputMenu
from ._shared.ca_value_float import CAValueFloat
from ._shared.ca_value_vector import CAValueVector
from .backbone_c import BackboneC
from .backbone_ca import BackboneCA
from .backbone_n import BackboneN
from .dihedral_angle import DihedralAngle
from .fallback_float import FallbackFloat
from .find_bonded_atom import FindBondedAtom
from .menu_residue_mask import MenuResidueMask


class DihedralPsi(AssetGeometryGroup):
    """
    Dihedral Psi

    Parameters
    ----------
    method : InputMenu | Literal["Read", "Compute"]
        Method

    Inputs
    ------
    i.method : MenuSocket
        Method

    Outputs
    -------
    o.psi : FloatSocket
        The calculated `Psi` angle for the residue, in the range of `(-pi, pi)`
    o.up : VectorSocket
        The perpendicular vector from the line of BC to the point A
    o.axis : VectorSocket
        The vector BC corresponding to the backbone vector around which the angle is calculated
    """

    _name = "Dihedral Psi"
    _asset_name = "Dihedral Psi"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.dihedral_psi"}

    class _Inputs(SocketAccessor):
        method: MenuSocket
        """Method"""

    class _Outputs(SocketAccessor):
        psi: FloatSocket
        """The calculated `Psi` angle for the residue, in the range of `(-pi, pi)`"""
        up: VectorSocket
        """The perpendicular vector from the line of BC to the point A"""
        axis: VectorSocket
        """The vector BC corresponding to the backbone vector around which the angle is calculated"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        method: InputMenu | Literal["Read", "Compute"] = "Compute",
    ):
        super().__init__(**{"Method": method})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        method = tree.inputs.menu("Method", expanded=True, optional_label=True)
        psi = tree.outputs.float(
            "Psi",
            description="The calculated `Psi` angle for the residue, in the range of `(-pi, pi)`",
            subtype="ANGLE",
        )
        up = tree.outputs.vector(
            "Up",
            description="The perpendicular vector from the line of BC to the point A",
        )
        axis = tree.outputs.vector(
            "Axis",
            description="The vector BC corresponding to the backbone vector around which the angle is calculated",
        )

        menu_switch = g.MenuSwitch.integer(method, {"Read": 0, "Compute": 1})
        group = DihedralAngle(
            a=BackboneN(
                method=g.IndexSwitch.menu(menu_switch.o.output, ("Read", "Compute"))
            ),
            b=BackboneCA(
                method=g.IndexSwitch.menu(menu_switch.o.output, ("Read", "Compute"))
            ),
            c=BackboneC(
                method=g.IndexSwitch.menu(menu_switch.o.output, ("Read", "Compute"))
            ),
            d=FindBondedAtom(
                index=MenuResidueMask(atom_name="C").o.index, distance=1
            ).o.position,
        )
        CAValueVector(vector=group.o.ba_bc) >> up
        CAValueVector(vector=group.o.bc) >> axis
        (
            CAValueFloat(
                value=FallbackFloat(name="Psi", fallback=group.o.angle)
            ).o.value
            * -1.0
            >> psi
        )

        method.default_value = "Compute"


ASSET = DihedralPsi

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
