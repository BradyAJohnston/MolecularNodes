# Node-group asset "PCA Dejitter" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import Default, InputGeometry, InputVector


class PCADejitter(AssetGeometryGroup):
    """
    Remove frame-to-frame jitter from a vector field during playback, by principal component analysis of its values over the last 4 frames

    Parameters
    ----------
    geometry : InputGeometry
        Points whose value is smoothed over time
    value : InputVector
        Vector field to remove jitter from, such as the position. When unconnected: The position from the context.

    Inputs
    ------
    i.geometry : GeometrySocket
        Points whose value is smoothed over time
    i.value : VectorSocket
        Vector field to remove jitter from, such as the position

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    o.value : VectorSocket
        The value with jitter removed, captured on the output geometry
    """

    _name = "PCA Dejitter"
    _asset_name = "PCA Dejitter"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Remove frame-to-frame jitter from a vector field during playback, by principal component analysis of its values over the last 4 frames"
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Points whose value is smoothed over time"""
        value: VectorSocket
        """Vector field to remove jitter from, such as the position"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        value: VectorSocket
        """The value with jitter removed, captured on the output geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        value: InputVector = Default.POSITION,
    ):
        super().__init__(**{"Geometry": geometry, "Value": value})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Points whose value is smoothed over time"
        )
        value = tree.inputs.vector(
            "Value",
            (0.0, 0.0, 0.0),
            description="Vector field to remove jitter from, such as the position",
            hide_value=True,
            default_input="POSITION",
        )
        geometry_1 = tree.outputs.geometry("Geometry")
        value_1 = tree.outputs.vector(
            "Value",
            description="The value with jitter removed, captured on the output geometry",
        )

        with g.Frame("History"):
            simulation_zone = g.SimulationZone()
            geometry_2 = simulation_zone.items.geometry(geometry, "Geometry")
            x0 = simulation_zone.items.vector(value, "X0")
            x1 = simulation_zone.items.vector(value, "X1")
            x2 = simulation_zone.items.vector(value, "X2")
            x3 = simulation_zone.items.vector(value, "X3")
            geometry >> geometry_2.next
            value >> x0.next
            (
                g.SampleIndex.point.vector(geometry_2.current, x0.current, g.Index())
                >> x1.next
            )
            (
                g.SampleIndex.point.vector(geometry_2.current, x1.current, g.Index())
                >> x2.next
            )
            (
                g.SampleIndex.point.vector(geometry_2.current, x2.current, g.Index())
                >> x3.next
            )
            _string = g.String(
                "Adapted from the PCA De-Jitter tutorial for Houdini by SideFX: https://www.sidefx.com/tutorials/pca-dejitter/ The simulation keeps the value of each point over the last 4 frames as X0 (this frame) to X3 (3 frames ago). Every step the history moves back one frame, read from the previous state by index, so the point count must not change. Scrubbing or jumping frames resets the history; it only builds up during playback."
            )
        with g.Frame("Covariance of the 4 frames"):
            combine_matrix = g.CombineMatrix(
                column_1_row_1=1.3333334,
                column_1_row_2=0.3333333,
                column_1_row_3=0.3333333,
                column_2_row_1=0.3333333,
                column_2_row_2=1.3333334,
                column_2_row_3=0.3333333,
                column_3_row_1=0.3333333,
                column_3_row_2=0.3333333,
                column_3_row_3=1.3333334,
            )
            _string_1 = g.String(
                "Each frame of the window is one PCA sample. Centred on the mean M, the samples D0..D3 sum to zero, so D1..D3 span them all and their 3x3 Gram matrix G (dot products summed over the points) holds the whole covariance. With K = I + 1/3 (all ones), the square root of the metric those three samples carry, H = K G K is symmetric and its SVD gives the principal components."
            )
            vector_math = (x0.result + x1.result + x2.result + x3.result) * 0.25
            vector_math_1 = x1.result - vector_math
            vector_math_2 = x2.result - vector_math
            vector_math_3 = x3.result - vector_math
            attribute_statistic = g.AttributeStatistic.point.float(
                geometry_2.result, attribute=vector_math_1.dot(vector_math_1)
            )
            attribute_statistic_1 = g.AttributeStatistic.point.float(
                geometry_2.result, attribute=vector_math_1.dot(vector_math_2)
            )
            attribute_statistic_2 = g.AttributeStatistic.point.float(
                geometry_2.result, attribute=vector_math_1.dot(vector_math_3)
            )
            attribute_statistic_3 = g.AttributeStatistic.point.float(
                geometry_2.result, attribute=vector_math_2.dot(vector_math_2)
            )
            attribute_statistic_4 = g.AttributeStatistic.point.float(
                geometry_2.result, attribute=vector_math_2.dot(vector_math_3)
            )
            attribute_statistic_5 = g.AttributeStatistic.point.float(
                geometry_2.result, attribute=vector_math_3.dot(vector_math_3)
            )
            combine_matrix_1 = g.CombineMatrix(
                column_1_row_1=attribute_statistic.o.sum,
                column_1_row_2=attribute_statistic_1.o.sum,
                column_1_row_3=attribute_statistic_2.o.sum,
                column_2_row_1=attribute_statistic_1.o.sum,
                column_2_row_2=attribute_statistic_3.o.sum,
                column_2_row_3=attribute_statistic_4.o.sum,
                column_3_row_1=attribute_statistic_2.o.sum,
                column_3_row_2=attribute_statistic_4.o.sum,
                column_3_row_3=attribute_statistic_5.o.sum,
            )
            multiply_matrices = g.MultiplyMatrices(
                matrix=g.MultiplyMatrices(
                    matrix=combine_matrix, matrix_001=combine_matrix_1
                ),
                matrix_001=combine_matrix,
            )
            matrix_svd = multiply_matrices.o.matrix.svd()
        with g.Frame("Fade components and rebuild"):
            _string_2 = g.String(
                "The current frame is projected onto the components and rebuilt with their weights faded from 1 (first) to 0.5 to 0 (last), as `molecularnodes.utils.pca_dejitter` does. The rebuilt value is M + B1 D1 + B2 D2 + B3 D3 with B = K U F U^T (-1/2, -1/2, -1/2). With every weight at 1 this gives back X0 exactly; at 0 it is the mean."
            )
            multiply_matrices_1 = g.MultiplyMatrices(
                matrix=matrix_svd.u,
                matrix_001=g.CombineMatrix(column_2_row_2=0.5, column_3_row_3=0.0),
            )
            multiply_matrices_2 = g.MultiplyMatrices(
                matrix=combine_matrix,
                matrix_001=g.MultiplyMatrices(
                    matrix=multiply_matrices_1, matrix_001=matrix_svd.u.transpose()
                ),
            )
            direction = g.TransformDirection(
                transform=multiply_matrices_2, direction=(-0.5, -0.5, -0.5)
            ).o.direction
            vector_math_4 = (
                vector_math
                + vector_math_1 * direction.x
                + vector_math_2 * direction.y
                + vector_math_3 * direction.z
            )
        capture = g.CaptureAttribute.point(geometry=geometry_2.result)
        value_2 = capture.items.vector(vector_math_4, "Value")

        capture.o.geometry >> geometry_1
        value_2.output >> value_1


ASSET = PCADejitter

ASSET_METADATA = {
    "description": "Remove frame-to-frame jitter from a vector field during playback, by principal component analysis of its values over the last 4 frames",
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
