# Node-group asset "Set Chi Angle" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .dihedral_chi_angle import DihedralChiAngle
from .peptide_chi import PeptideChi
from .set_ures_id import SetUResID


class SetChiAngle(AssetGeometryGroup):
    """
    Set Chi Angle

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    selection : InputBoolean
        Selection
    x1 : InputFloat
        X1
    x2 : InputFloat
        X2
    x3 : InputFloat
        X3
    x4 : InputFloat
        X4
    x5 : InputFloat
        X5

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.selection : BooleanSocket
        Selection
    i.x1 : FloatSocket
        X1
    i.x2 : FloatSocket
        X2
    i.x3 : FloatSocket
        X3
    i.x4 : FloatSocket
        X4
    i.x5 : FloatSocket
        X5

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Set Chi Angle"
    _asset_name = "Set Chi Angle"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        selection: BooleanSocket
        """Selection"""
        x1: FloatSocket
        """X1"""
        x2: FloatSocket
        """X2"""
        x3: FloatSocket
        """X3"""
        x4: FloatSocket
        """X4"""
        x5: FloatSocket
        """X5"""

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
        x1: InputFloat = 0.0,
        x2: InputFloat = 0.0,
        x3: InputFloat = 0.0,
        x4: InputFloat = 0.0,
        x5: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Selection": selection,
                "X1": x1,
                "X2": x2,
                "X3": x3,
                "X4": x4,
                "X5": x5,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        x1 = tree.inputs.float(
            "X1", 0.0, min_value=-10_000.0, max_value=10_000.0, subtype="ANGLE"
        )
        x2 = tree.inputs.float(
            "X2", 0.0, min_value=-10_000.0, max_value=10_000.0, subtype="ANGLE"
        )
        x3 = tree.inputs.float(
            "X3", 0.0, min_value=-10_000.0, max_value=10_000.0, subtype="ANGLE"
        )
        x4 = tree.inputs.float(
            "X4", 0.0, min_value=-10_000.0, max_value=10_000.0, subtype="ANGLE"
        )
        x5 = tree.inputs.float(
            "X5", 0.0, min_value=-10_000.0, max_value=10_000.0, subtype="ANGLE"
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        math_1 = DihedralChiAngle().o.angle * -1.0
        group = PeptideChi(
            selection=selection,
            x1=x1 + math_1,
            x2=x2 + math_1,
            x3=x3 + math_1,
            x4=x4 + math_1,
            x5=x5 + math_1,
        )
        SetUResID(geometry=geometry) >> g.SetPosition(position=group) >> geometry_1


ASSET = SetChiAngle

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
