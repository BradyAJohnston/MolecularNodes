# Node-group asset "Mesh Contours" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    FloatSocket,
    GeometrySocket,
    MaterialSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputFloat, InputGeometry, InputMaterial
from .set_color import SetColor


class MeshContours(AssetGeometryGroup):
    """
    Mesh Contours

    Parameters
    ----------
    mesh : InputGeometry
        ISO surface geometry output
    value : InputFloat
        Value
    color : InputColor
        Color to apply to the selected atoms
    material : InputMaterial
        Material

    Inputs
    ------
    i.mesh : GeometrySocket
        ISO surface geometry output
    i.value : FloatSocket
        Value
    i.color : ColorSocket
        Color to apply to the selected atoms
    i.material : MaterialSocket
        Material

    Outputs
    -------
    o.mesh : GeometrySocket
        Mesh
    """

    _name = "Mesh Contours"
    _asset_name = "Mesh Contours"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        mesh: GeometrySocket
        """ISO surface geometry output"""
        value: FloatSocket
        """Value"""
        color: ColorSocket
        """Color to apply to the selected atoms"""
        material: MaterialSocket
        """Material"""

    class _Outputs(SocketAccessor):
        mesh: GeometrySocket
        """Mesh"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        mesh: InputGeometry = None,
        value: InputFloat = 1.0,
        color: InputColor = None,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{"Mesh": mesh, "Value": value, "Color": color, "Material": material}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        mesh = tree.inputs.geometry("Mesh", description="ISO surface geometry output")
        value = tree.inputs.float("Value", 1.0, min_value=-10_000.0, max_value=10_000.0)
        color = tree.inputs.color(
            "Color",
            (0.0, 0.0, 0.0, 1.0),
            description="Color to apply to the selected atoms",
        )
        material = tree.inputs.material("Material", optional_label=True)
        mesh_1 = tree.outputs.geometry("Mesh")

        math_1 = value / 1000.0
        set_material = (
            mesh
            >> g.SplitEdges()
            >> g.MeshToCurve()
            >> g.CurveToMesh(profile_curve=g.Quadrilateral(width=math_1, height=math_1))
            >> g.SetMaterial(material=material)
        )
        SetColor(atoms=set_material, color=color) >> mesh_1


ASSET = MeshContours

ASSET_METADATA = {
    "description": "Turn the edges of a mesh into mesh tubes. Useful for styling the contours of a surface mesh as wires.",
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
