# Node-group asset '3 Index Angle' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .vector_angle import VectorAngle


class Group3IndexAngle(AssetGeometryGroup):
    """
    3 Index Angle

    Parameters
    ----------
    index_a : InputInteger
        First of the points for the angle calculation
    index_b : InputInteger
        The middle point for the angle calculation
    index_c : InputInteger
        Last of the points for the angle calculation

    Inputs
    ------
    i.index_a : IntegerSocket
        First of the points for the angle calculation
    i.index_b : IntegerSocket
        The middle point for the angle calculation
    i.index_c : IntegerSocket
        Last of the points for the angle calculation

    Outputs
    -------
    o.angle : FloatSocket
        Angle between the points around Index B in radians
    """

    _name = "3 Index Angle"
    _asset_name = "3 Index Angle"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.3_point_angle"}

    class _Inputs(SocketAccessor):
        index_a: IntegerSocket
        """First of the points for the angle calculation"""
        index_b: IntegerSocket
        """The middle point for the angle calculation"""
        index_c: IntegerSocket
        """Last of the points for the angle calculation"""

    class _Outputs(SocketAccessor):
        angle: FloatSocket
        """Angle between the points around Index B in radians"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index_a: InputInteger = 0,
        index_b: InputInteger = 1,
        index_c: InputInteger = 2,
    ):
        super().__init__(**{"Index A": index_a, "Index B": index_b, "Index C": index_c})

    def _build_group(self, tree):
        index_a = tree.inputs.integer(
            "Index A",
            0,
            description="First of the points for the angle calculation",
            min_value=0,
        )
        index_b = tree.inputs.integer(
            "Index B",
            1,
            description="The middle point for the angle calculation",
            min_value=0,
        )
        index_c = tree.inputs.integer(
            "Index C",
            2,
            description="Last of the points for the angle calculation",
            min_value=0,
        )
        angle = tree.outputs.float(
            "Angle",
            description="Angle between the points around Index B in radians",
            subtype="ANGLE",
        )

        position = g.Position()
        evaluate_at_index = position.o.position.point.at(index_b)
        (
            VectorAngle(
                a=position.o.position.point.at(index_a) - evaluate_at_index,
                b=position.o.position.point.at(index_c) - evaluate_at_index,
            )
            >> angle
        )


ASSET = Group3IndexAngle

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
