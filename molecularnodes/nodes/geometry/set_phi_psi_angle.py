# Node-group asset "Set Phi Psi Angle" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .dihedral_phi import DihedralPhi
from .dihedral_psi import DihedralPsi
from .peptide_dihedral import PeptideDihedral


class SetPhiPsiAngle(AssetGeometryGroup):
    """
    Set Phi Psi Angle

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    selection : InputBoolean
        Selection
    phi : InputFloat
        Phi
    psi : InputFloat
        Psi

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.selection : BooleanSocket
        Selection
    i.phi : FloatSocket
        Phi
    i.psi : FloatSocket
        Psi

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Set Phi Psi Angle"
    _asset_name = "Set Phi Psi Angle"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        selection: BooleanSocket
        """Selection"""
        phi: FloatSocket
        """Phi"""
        psi: FloatSocket
        """Psi"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        selection: InputBoolean = True,
        phi: InputFloat = 0.0,
        psi: InputFloat = 0.0,
    ):
        super().__init__(
            **{"Geometry": geometry, "Selection": selection, "Phi": phi, "Psi": psi}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        phi = tree.inputs.float(
            "Phi", 0.0, min_value=-10_000.0, max_value=10_000.0, subtype="ANGLE"
        )
        psi = tree.inputs.float(
            "Psi", 0.0, min_value=-10_000.0, max_value=10_000.0, subtype="ANGLE"
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        group = PeptideDihedral(
            selection=selection,
            phi=DihedralPhi().o.phi.mul_add(-1.0, phi),
            psi=DihedralPsi().o.psi.mul_add(-1.0, psi),
        )
        geometry >> g.SetPosition(position=group) >> geometry_1


ASSET = SetPhiPsiAngle

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
