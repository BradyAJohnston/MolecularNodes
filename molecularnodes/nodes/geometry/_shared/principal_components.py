# Node group 'Principal Components' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    IntegerSocket,
    RotationSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputInteger, InputVector


class PrincipalComponents(CustomGeometryGroup):
    """
    Principal Component Analysis (PCA) of a vector field

    Parameters
    ----------
    position : InputVector
        Position
    group_id : InputInteger
        An index used to group values together for multiple separate operations

    Inputs
    ------
    i.position : VectorSocket
        Position
    i.group_id : IntegerSocket
        An index used to group values together for multiple separate operations

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

    _name = "Principal Components"
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Principal Component Analysis (PCA) of a vector field",
        "default_group_node_width": 200,
    }

    class _Inputs(SocketAccessor):
        position: VectorSocket
        """Position"""
        group_id: IntegerSocket
        """An index used to group values together for multiple separate operations"""

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
        position: InputVector = None,
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Position": position, "Group ID": group_id})

    def _build_group(self, tree):
        position = tree.inputs.vector(
            "Position", (0.0, 0.0, 0.0), subtype="TRANSLATION", default_input="POSITION"
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="An index used to group values together for multiple separate operations",
            hide_value=True,
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

        with g.Frame("Centroid"):
            field_average = position.point.mean(group_id)
        with g.Frame("Covariance matrix"):
            vector_math = position - field_average
            mean = (vector_math * vector_math.x).point.mean(group_id)
            mean_1 = (vector_math * vector_math.y).point.mean(group_id)
            mean_2 = (vector_math * vector_math.z).point.mean(group_id)
            combine_matrix = g.CombineMatrix(
                column_1_row_1=mean.x,
                column_1_row_2=mean.y,
                column_1_row_3=mean.z,
                column_2_row_1=mean_1.x,
                column_2_row_2=mean_1.y,
                column_2_row_3=mean_1.z,
                column_3_row_1=mean_2.x,
                column_3_row_2=mean_2.y,
                column_3_row_3=mean_2.z,
            )
        with g.Frame("SVD of Covariance to find principal components"):
            matrix_svd = combine_matrix.o.matrix.svd()
            matrix_determinant = matrix_svd.u.determinant()
        separate_matrix = g.SeparateMatrix(matrix=matrix_svd.u)
        combine_xyz = g.CombineXYZ(
            x=separate_matrix.o.column_1_row_1,
            y=separate_matrix.o.column_1_row_2,
            z=separate_matrix.o.column_1_row_3,
        )
        combine_xyz_1 = g.CombineXYZ(
            x=separate_matrix.o.column_2_row_1,
            y=separate_matrix.o.column_2_row_2,
            z=separate_matrix.o.column_2_row_3,
        )
        combine_xyz_2 = g.CombineXYZ(
            x=separate_matrix.o.column_3_row_1,
            y=separate_matrix.o.column_3_row_2,
            z=separate_matrix.o.column_3_row_3,
        )
        axes_to_rotation = g.AxesToRotation(
            primary_axis=combine_xyz,
            secondary_axis=combine_xyz_2,
            primary="X",
            secondary="Z",
        )
        combine_xyz_1.o.vector * matrix_determinant.sign() >> intermediate_axis

        field_average >> group_center
        axes_to_rotation >> rotation
        matrix_svd.s >> principal_components
        combine_xyz >> longest_axis
        combine_xyz_2 >> shortest_axis
