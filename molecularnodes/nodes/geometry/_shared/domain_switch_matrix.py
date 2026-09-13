# Node group "Domain Switch Matrix" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    MatrixSocket,
    MenuSocket,
    SocketAccessor,
)
from nodebpy.types import InputMatrix, InputMenu


class DomainSwitchMatrix(CustomGeometryGroup):
    """
    Domain Switch Matrix

    Parameters
    ----------
    menu : InputMenu | Literal["Point", "Edge", "Face", "Face Corner", "Spline", "Instance"]
        Menu
    point : InputMatrix
        Point
    edge : InputMatrix
        Edge
    face : InputMatrix
        Face
    face_corner : InputMatrix
        Face Corner
    spline : InputMatrix
        Spline
    instance : InputMatrix
        Instance

    Inputs
    ------
    i.menu : MenuSocket
        Menu
    i.point : MatrixSocket
        Point
    i.edge : MatrixSocket
        Edge
    i.face : MatrixSocket
        Face
    i.face_corner : MatrixSocket
        Face Corner
    i.spline : MatrixSocket
        Spline
    i.instance : MatrixSocket
        Instance

    Outputs
    -------
    o.output : MatrixSocket
        Output
    """

    _name = "Domain Switch Matrix"
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.domain_switch_matrix"}

    class _Inputs(SocketAccessor):
        menu: MenuSocket
        """Menu"""
        point: MatrixSocket
        """Point"""
        edge: MatrixSocket
        """Edge"""
        face: MatrixSocket
        """Face"""
        face_corner: MatrixSocket
        """Face Corner"""
        spline: MatrixSocket
        """Spline"""
        instance: MatrixSocket
        """Instance"""

    class _Outputs(SocketAccessor):
        output: MatrixSocket
        """Output"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        menu: InputMenu
        | Literal[
            "Point", "Edge", "Face", "Face Corner", "Spline", "Instance"
        ] = "Point",
        point: InputMatrix = None,
        edge: InputMatrix = None,
        face: InputMatrix = None,
        face_corner: InputMatrix = None,
        spline: InputMatrix = None,
        instance: InputMatrix = None,
    ):
        super().__init__(
            **{
                "Menu": menu,
                "Point": point,
                "Edge": edge,
                "Face": face,
                "Face Corner": face_corner,
                "Spline": spline,
                "Instance": instance,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        menu = tree.inputs.menu("Menu", optional_label=True)
        point = tree.inputs.matrix("Point")
        edge = tree.inputs.matrix("Edge")
        face = tree.inputs.matrix("Face")
        face_corner = tree.inputs.matrix("Face Corner")
        spline = tree.inputs.matrix("Spline")
        instance = tree.inputs.matrix("Instance")
        output = tree.outputs.matrix("Output")

        (
            g.MenuSwitch.matrix(
                menu,
                {
                    "Point": point,
                    "Edge": edge,
                    "Face": face,
                    "Face Corner": face_corner,
                    "Spline": spline,
                    "Instance": instance,
                },
            )
            >> output
        )

        menu.default_value = "Point"
