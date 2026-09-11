# Node-group asset "Curve Offset Dihedral" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputInteger, InputVector
from .dihedral_angle import DihedralAngle
from .offset_vector import OffsetVector


class CurveOffsetDihedral(AssetGeometryGroup):
    """
    Offset from the current point a number of points, then use their `Position` and `Normal` to calculat a dihdral angle between them

    Parameters
    ----------
    position : InputVector
        The vector to use as the B & C components for `Dihedral Angle` calculation
    normal : InputVector
        The normal that will be added the `Position` to create the A & D components of the dihedral calcaulation
    index : InputInteger
        The index of the current point to calculate from
    offset : InputInteger
        The number of points to offset before calculating the angle

    Inputs
    ------
    i.position : VectorSocket
        The vector to use as the B & C components for `Dihedral Angle` calculation
    i.normal : VectorSocket
        The normal that will be added the `Position` to create the A & D components of the dihedral calcaulation
    i.index : IntegerSocket
        The index of the current point to calculate from
    i.offset : IntegerSocket
        The number of points to offset before calculating the angle

    Outputs
    -------
    o.angle : FloatSocket
        The calculated angle in radians
    """

    _name = "Curve Offset Dihedral"
    _asset_name = "Curve Offset Dihedral"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Offset from the current point a number of points, then use their `Position` and `Normal` to calculat a dihdral angle between them",
        "node_tool_idname": "geometry.curve_offset_dihedral",
    }

    class _Inputs(SocketAccessor):
        position: VectorSocket
        """The vector to use as the B & C components for `Dihedral Angle` calculation"""
        normal: VectorSocket
        """The normal that will be added the `Position` to create the A & D components of the dihedral calcaulation"""
        index: IntegerSocket
        """The index of the current point to calculate from"""
        offset: IntegerSocket
        """The number of points to offset before calculating the angle"""

    class _Outputs(SocketAccessor):
        angle: FloatSocket
        """The calculated angle in radians"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        position: InputVector = None,
        normal: InputVector = None,
        index: InputInteger = 0,
        offset: InputInteger = 0,
    ):
        super().__init__(
            **{"Position": position, "Normal": normal, "Index": index, "Offset": offset}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        position = tree.inputs.vector(
            "Position",
            (0.0, 0.0, 0.0),
            description="The vector to use as the B & C components for `Dihedral Angle` calculation",
            default_input="POSITION",
        )
        normal = tree.inputs.vector(
            "Normal",
            (0.0, 0.0, 0.0),
            description="The normal that will be added the `Position` to create the A & D components of the dihedral calcaulation",
            default_input="NORMAL",
        )
        index = tree.inputs.integer(
            "Index",
            0,
            description="The index of the current point to calculate from",
            min_value=-2147483647,
            default_input="INDEX",
        )
        offset = tree.inputs.integer(
            "Offset",
            0,
            description="The number of points to offset before calculating the angle",
        )
        angle = tree.outputs.float(
            "Angle", description="The calculated angle in radians", subtype="ANGLE"
        )

        group = OffsetVector(vector=position, index=index, offset=offset)
        (
            DihedralAngle(
                a=OffsetVector(vector=normal, index=index, offset=offset).o.value
                + group,
                b=group,
                c=position,
                d=position + normal,
            )
            >> angle
        )


ASSET = CurveOffsetDihedral

ASSET_METADATA = {
    "description": "Offset from the current point a number of points, then use their `Position` and `Normal` to calculat a dihdral angle between them",
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
