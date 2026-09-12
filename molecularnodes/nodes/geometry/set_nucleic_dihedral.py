# Node-group asset "Set Nucleic Dihedral" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .dihedral_nucleic_angle import DihedralNucleicAngle
from .nucleic_dihedral import NucleicDihedral


class SetNucleicDihedral(AssetGeometryGroup):
    """
    Set Nucleic Dihedral

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    selection : InputBoolean
        The resulting selection must overlap with this input selection
    alpha : InputFloat
        Alpha
    beta : InputFloat
        Beta
    gamma : InputFloat
        Gamma
    epsilon : InputFloat
        Epsilon
    zeta : InputFloat
        Amount to rotate around the axis

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.selection : BooleanSocket
        The resulting selection must overlap with this input selection
    i.alpha : FloatSocket
        Alpha
    i.beta : FloatSocket
        Beta
    i.gamma : FloatSocket
        Gamma
    i.epsilon : FloatSocket
        Epsilon
    i.zeta : FloatSocket
        Amount to rotate around the axis

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Set Nucleic Dihedral"
    _asset_name = "Set Nucleic Dihedral"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        selection: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        alpha: FloatSocket
        """Alpha"""
        beta: FloatSocket
        """Beta"""
        gamma: FloatSocket
        """Gamma"""
        epsilon: FloatSocket
        """Epsilon"""
        zeta: FloatSocket
        """Amount to rotate around the axis"""

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
        alpha: InputFloat = 0.0,
        beta: InputFloat = 0.0,
        gamma: InputFloat = 0.0,
        epsilon: InputFloat = 0.0,
        zeta: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Selection": selection,
                "Alpha": alpha,
                "Beta": beta,
                "Gamma": gamma,
                "Epsilon": epsilon,
                "Zeta": zeta,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="The resulting selection must overlap with this input selection",
            hide_value=True,
        )
        alpha = tree.inputs.float("Alpha", 0.0, subtype="ANGLE")
        beta = tree.inputs.float("Beta", 0.0, subtype="ANGLE")
        gamma = tree.inputs.float("Gamma", 0.0, subtype="ANGLE")
        epsilon = tree.inputs.float("Epsilon", 0.0, subtype="ANGLE")
        zeta = tree.inputs.float(
            "Zeta", 0.0, description="Amount to rotate around the axis", subtype="ANGLE"
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        math_1 = DihedralNucleicAngle().o.angle * -1.0
        group = NucleicDihedral(
            selection=selection,
            alpha=alpha + math_1,
            beta=beta + math_1,
            gamma=gamma + math_1,
            epsilon=epsilon + math_1,
            zeta=zeta + math_1,
        )
        geometry >> g.SetPosition(position=group) >> geometry_1


ASSET = SetNucleicDihedral

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
