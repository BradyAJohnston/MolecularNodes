# Node-group asset "Dihedral Angle" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputVector
from .vector_angle import VectorAngle


class DihedralAngle(AssetGeometryGroup):
    """
    Dihedral Angle

    Parameters
    ----------
    a : InputVector
        First vector for the calculation, which draws a line to B
    b : InputVector
        Second vector for the calculation, which receives a line from A and draws a line to C
    c : InputVector
        Third vector for the calculation, which receives a line from B and draws a line to D
    d : InputVector
        Last vector for the calculation, which is the end point of the line from D

    Inputs
    ------
    i.a : VectorSocket
        First vector for the calculation, which draws a line to B
    i.b : VectorSocket
        Second vector for the calculation, which receives a line from A and draws a line to C
    i.c : VectorSocket
        Third vector for the calculation, which receives a line from B and draws a line to D
    i.d : VectorSocket
        Last vector for the calculation, which is the end point of the line from D

    Outputs
    -------
    o.angle : FloatSocket
        The angle between the vectors AB and CD, when made perpendicular to BC.
    o.ba_bc : VectorSocket
        The vector BA when made perpendicular to the axis BC
    o.cd_bc : VectorSocket
        The Vector CD when makde perpendicular to the axis BC
    o.bc : VectorSocket
        The axis vector BC
    """

    _name = "Dihedral Angle"
    _asset_name = "Dihedral Angle"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "VECTOR"
    _tree_properties = {"node_tool_idname": "geometry.dihedral_angle"}

    class _Inputs(SocketAccessor):
        a: VectorSocket
        """First vector for the calculation, which draws a line to B"""
        b: VectorSocket
        """Second vector for the calculation, which receives a line from A and draws a line to C"""
        c: VectorSocket
        """Third vector for the calculation, which receives a line from B and draws a line to D"""
        d: VectorSocket
        """Last vector for the calculation, which is the end point of the line from D"""

    class _Outputs(SocketAccessor):
        angle: FloatSocket
        """The angle between the vectors AB and CD, when made perpendicular to BC."""
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
        a: InputVector = None,
        b: InputVector = None,
        c: InputVector = None,
        d: InputVector = None,
    ):
        super().__init__(**{"A": a, "B": b, "C": c, "D": d})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        a = tree.inputs.vector(
            "A",
            (0.0, 0.0, 0.0),
            description="First vector for the calculation, which draws a line to B",
        )
        b = tree.inputs.vector(
            "B",
            (0.0, 0.0, 0.0),
            description="Second vector for the calculation, which receives a line from A and draws a line to C",
        )
        c_ = tree.inputs.vector(
            "C",
            (0.0, 0.0, 0.0),
            description="Third vector for the calculation, which receives a line from B and draws a line to D",
        )
        d = tree.inputs.vector(
            "D",
            (0.0, 0.0, 0.0),
            description="Last vector for the calculation, which is the end point of the line from D",
        )
        angle = tree.outputs.float(
            "Angle",
            description="The angle between the vectors AB and CD, when made perpendicular to BC.",
            subtype="ANGLE",
        )
        ba_bc = tree.outputs.vector(
            "BA⟂(BC)",
            description="The vector BA when made perpendicular to  the axis BC",
        )
        cd_bc = tree.outputs.vector(
            "CD⟂(BC)",
            description="The Vector CD when makde perpendicular to the axis BC",
        )
        bc = tree.outputs.vector("BC", description="The axis vector BC")

        vector_math = d - c_
        vector_math_1 = a - b
        vector_math_2 = c_ - b
        vector_math_3 = vector_math_1 - vector_math_1.project(vector_math_2)
        vector_math_4 = vector_math - vector_math.project(vector_math_2)
        (
            VectorAngle(a=vector_math_3, b=vector_math_4).o.angle
            * vector_math_3.cross(vector_math_4).dot(vector_math_2).sign()
            >> angle
        )

        vector_math_3 >> ba_bc
        vector_math_4 >> cd_bc
        vector_math_2 >> bc


ASSET = DihedralAngle

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
