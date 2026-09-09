# Node group '.OKLab Matrices' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, MatrixSocket, SocketAccessor


class OKLabMatrices(CustomGeometryGroup):
    """
    .OKLab Matrices

    Outputs
    -------
    o.m1 : MatrixSocket
        M1
    o.m2 : MatrixSocket
        M2
    """

    _name = ".OKLab Matrices"
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        m1: MatrixSocket
        """M1"""
        m2: MatrixSocket
        """M2"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        m1 = tree.outputs.matrix("M1")
        m2 = tree.outputs.matrix("M2")

        combine_matrix = g.CombineMatrix(
            column_1_row_1=0.818933,
            column_1_row_2=0.032985,
            column_1_row_3=0.0482,
            column_2_row_1=0.3618667,
            column_2_row_2=0.9293119,
            column_2_row_3=0.264366,
            column_3_row_1=-0.1288597,
            column_3_row_2=0.03614564,
            column_3_row_3=0.633852,
        )
        combine_matrix_1 = g.CombineMatrix(
            column_1_row_1=0.21045426,
            column_1_row_2=1.9779985,
            column_1_row_3=0.025904037,
            column_2_row_1=0.7936178,
            column_2_row_2=-2.42859,
            column_2_row_3=0.7827718,
            column_3_row_1=-0.004072047,
            column_3_row_2=0.4505937,
            column_3_row_3=-0.8086758,
        )

        combine_matrix >> m1
        combine_matrix_1 >> m2
