# Node group "Geometry Principal Components" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    GeometrySocket,
    RotationSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputGeometry, InputVector
from .principal_components import PrincipalComponents


class GeometryPrincipalComponents(CustomGeometryGroup):
    """
    Principal Component Analysis (PCA) of a the positions of points

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to evaluate the given fields and store the resulting attributes on. All geometry types except volumes are supported
    position : InputVector
        Position

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to evaluate the given fields and store the resulting attributes on. All geometry types except volumes are supported
    i.position : VectorSocket
        Position

    Outputs
    -------
    o.group_center : VectorSocket
        Group Center
    o.rotation : RotationSocket
        Rotation that defines the principal component basis
    o.principal_components : VectorSocket
        Variance of the data along each principal axis
    o.longest_axis : VectorSocket
        Axis along the most variance
    o.intermediate_axis : VectorSocket
        Axis completing a right handed orthogonal basis with the longest and shortest axis
    o.shortest_axis : VectorSocket
        Axis along the least variance
    """

    _name = "Geometry Principal Components"
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Principal Component Analysis (PCA) of a the positions of points",
        "default_group_node_width": 200,
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to evaluate the given fields and store the resulting attributes on. All geometry types except volumes are supported"""
        position: VectorSocket
        """Position"""

    class _Outputs(SocketAccessor):
        group_center: VectorSocket
        """Group Center"""
        rotation: RotationSocket
        """Rotation that defines the principal component basis"""
        principal_components: VectorSocket
        """Variance of the data along each principal axis"""
        longest_axis: VectorSocket
        """Axis along the most variance"""
        intermediate_axis: VectorSocket
        """Axis completing a right handed orthogonal basis with the longest and shortest axis"""
        shortest_axis: VectorSocket
        """Axis along the least variance"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        position: InputVector = None,
    ):
        super().__init__(**{"Geometry": geometry, "Position": position})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry",
            description="Geometry to evaluate the given fields and store the resulting attributes on. All geometry types except volumes are supported",
        )
        position = tree.inputs.vector(
            "Position", (0.0, 0.0, 0.0), subtype="TRANSLATION", default_input="POSITION"
        )
        group_center = tree.outputs.vector("Group Center")
        rotation = tree.outputs.rotation(
            "Rotation",
            description="Rotation that defines the principal component basis",
        )
        principal_components = tree.outputs.vector(
            "Principal Components",
            description="Variance of the data along each principal axis",
        )
        with tree.outputs.panel("Principal Axes"):
            longest_axis = tree.outputs.vector(
                "Longest Axis", description="Axis along the most variance"
            )
            intermediate_axis = tree.outputs.vector(
                "Intermediate Axis",
                description="Axis completing a right handed orthogonal basis with the longest and shortest axis",
            )
            shortest_axis = tree.outputs.vector(
                "Shortest Axis", description="Axis along the least variance"
            )

        capture = g.CaptureAttribute.point(geometry=geometry)
        position_1 = capture.items.vector("Position", position)
        group = PrincipalComponents(position=position_1.output)
        sample_index = g.SampleIndex(
            geometry=capture.o.geometry,
            value=group.o.group_center,
            data_type="FLOAT_VECTOR",
        )
        sample_index_1 = g.SampleIndex(
            geometry=capture.o.geometry, value=group.o.rotation, data_type="QUATERNION"
        )
        sample_index_2 = g.SampleIndex(
            geometry=capture.o.geometry,
            value=group.o.principal_components,
            data_type="FLOAT_VECTOR",
        )
        sample_index_3 = g.SampleIndex(
            geometry=capture.o.geometry,
            value=group.o.longest_axis,
            data_type="FLOAT_VECTOR",
        )
        sample_index_4 = g.SampleIndex(
            geometry=capture.o.geometry,
            value=group.o.intermediate_axis,
            data_type="FLOAT_VECTOR",
        )
        sample_index_5 = g.SampleIndex(
            geometry=capture.o.geometry,
            value=group.o.shortest_axis,
            data_type="FLOAT_VECTOR",
        )

        sample_index >> group_center
        sample_index_1 >> rotation
        sample_index_2 >> principal_components
        sample_index_3 >> longest_axis
        sample_index_4 >> intermediate_axis
        sample_index_5 >> shortest_axis
