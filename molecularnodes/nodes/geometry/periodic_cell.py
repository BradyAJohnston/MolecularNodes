# Node-group asset "Periodic Cell" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputMenu
from ._shared.periodic_nearest_image import PeriodicNearestImage
from .periodic_box import PeriodicBox


class PeriodicCell(AssetGeometryGroup):
    """
    Mesh and edges of the periodic unit cell, the box spanned by the box vectors or the compact cell around its centre

    Parameters
    ----------
    wrap : InputMenu | Literal["Box", "Compact"]
        The box spanned by the box vectors, or the compact cell around the box centre (GROMACS `-ur compact`), matching the Wrap of Periodic Offset
    update : InputBoolean
        Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node. As a modifier input, reads the "True" attribute by default.
    a : InputFloat
        a
    b : InputFloat
        b
    c : InputFloat
        c
    alpha : InputFloat
        alpha
    beta : InputFloat
        beta
    gamma : InputFloat
        gamma

    Inputs
    ------
    i.wrap : MenuSocket
        The box spanned by the box vectors, or the compact cell around the box centre (GROMACS `-ur compact`), matching the Wrap of Periodic Offset
    i.update : BooleanSocket
        Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node.
    i.a : FloatSocket
        a
    i.b : FloatSocket
        b
    i.c : FloatSocket
        c
    i.alpha : FloatSocket
        alpha
    i.beta : FloatSocket
        beta
    i.gamma : FloatSocket
        gamma

    Outputs
    -------
    o.mesh : GeometrySocket
        The closed cell, empty without a box
    o.edges : GeometrySocket
        The edges of the cell as curves, without the diagonals of its triangulated faces
    """

    _name = "Periodic Cell"
    _asset_name = "Periodic Cell"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Mesh and edges of the periodic unit cell, the box spanned by the box vectors or the compact cell around its centre"
    }

    class _Inputs(SocketAccessor):
        wrap: MenuSocket
        """The box spanned by the box vectors, or the compact cell around the box centre (GROMACS `-ur compact`), matching the Wrap of Periodic Offset"""
        update: BooleanSocket
        """Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node."""
        a: FloatSocket
        b: FloatSocket
        c: FloatSocket
        alpha: FloatSocket
        beta: FloatSocket
        gamma: FloatSocket

    class _Outputs(SocketAccessor):
        mesh: GeometrySocket
        """The closed cell, empty without a box"""
        edges: GeometrySocket
        """The edges of the cell as curves, without the diagonals of its triangulated faces"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        wrap: InputMenu | Literal["Box", "Compact"] = "Box",
        update: InputBoolean = True,
        a: InputFloat = 0.0,
        b: InputFloat = 0.0,
        c: InputFloat = 0.0,
        alpha: InputFloat = 0.0,
        beta: InputFloat = 0.0,
        gamma: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Wrap": wrap,
                "Update": update,
                "a": a,
                "b": b,
                "c": c,
                "alpha": alpha,
                "beta": beta,
                "gamma": gamma,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        wrap = tree.inputs.menu(
            "Wrap",
            "Box",
            description="The box spanned by the box vectors, or the compact cell around the box centre (GROMACS `-ur compact`), matching the Wrap of Periodic Offset",
            expanded=True,
            optional_label=True,
        )
        with tree.inputs.panel("Periodic Box", default_closed=True):
            _update = tree.inputs.boolean(
                "Update",
                True,
                description="Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node.",
                default_attribute="True",
            )
            with tree.inputs.panel("Lengths"):
                a = tree.inputs.float("a", 0.0)
                b = tree.inputs.float("b", 0.0)
                c_ = tree.inputs.float("c", 0.0)
            with tree.inputs.panel("Angles"):
                alpha = tree.inputs.float("alpha", 0.0)
                beta = tree.inputs.float("beta", 0.0)
                gamma = tree.inputs.float("gamma", 0.0)
        mesh = tree.outputs.geometry(
            "Mesh", description="The closed cell, empty without a box"
        )
        edges = tree.outputs.geometry(
            "Edges",
            description="The edges of the cell as curves, without the diagonals of its triangulated faces",
        )

        periodic_box = PeriodicBox(a=a, b=b, c=c_, alpha=alpha, beta=beta, gamma=gamma)
        with g.Frame("Box"):
            index = g.Index()
            _string = g.String("The 8 corners of the box spanned by the box vectors.")
            transform_point = g.CombineXYZ(
                x=index.o.index % 2, y=index.o.index // 2 % 2, z=index.o.index // 4
            ).o.vector.transform(periodic_box.o.matrix)
            points = g.Points(position=transform_point, count=8)
        with g.Frame("Compact"):
            index_1 = g.Index()
            integer_math = index_1.o.index % 26
            integer_math_1 = integer_math + integer_math // 13
            integer_math_2 = index_1.o.index // 26 % 26
            integer_math_3 = integer_math_2 + integer_math_2 // 13
            integer_math_4 = index_1.o.index // 676
            integer_math_5 = integer_math_4 + integer_math_4 // 13
            _string_1 = g.String(
                "Each corner of the compact cell lies on three of the planes half way to the 26 neighbouring lattice points. Solve for the point on every three planes, then keep those with no nearer lattice point than the centre, which are the corners."
            )
            vector_math = periodic_box.o.a * (
                integer_math_1 % 3 - 1
            ) + periodic_box.o.b * (integer_math_1 // 3 % 3 - 1)
            vector_math_1 = vector_math + periodic_box.o.c * (integer_math_1 // 9 - 1)
            vector_math_2 = periodic_box.o.a * (
                integer_math_3 % 3 - 1
            ) + periodic_box.o.b * (integer_math_3 // 3 % 3 - 1)
            vector_math_3 = vector_math_2 + periodic_box.o.c * (integer_math_3 // 9 - 1)
            vector_math_4 = periodic_box.o.a * (
                integer_math_5 % 3 - 1
            ) + periodic_box.o.b * (integer_math_5 // 3 % 3 - 1)
            vector_math_5 = vector_math_4 + periodic_box.o.c * (integer_math_5 // 9 - 1)
            combine_matrix = g.CombineMatrix(
                column_1_row_1=vector_math_1.x,
                column_1_row_2=vector_math_3.x,
                column_1_row_3=vector_math_5.x,
                column_2_row_1=vector_math_1.y,
                column_2_row_2=vector_math_3.y,
                column_2_row_3=vector_math_5.y,
                column_3_row_1=vector_math_1.z,
                column_3_row_2=vector_math_3.z,
                column_3_row_3=vector_math_5.z,
            )
            invert_matrix = g.InvertMatrix(combine_matrix)
            combine_xyz = g.CombineXYZ(
                x=vector_math_1.dot(vector_math_1) * 0.5,
                y=vector_math_3.dot(vector_math_3) * 0.5,
                z=vector_math_5.dot(vector_math_5) * 0.5,
            )
            transform_point_1 = combine_xyz.o.vector.transform(invert_matrix.o.matrix)
            periodic_nearest_image = PeriodicNearestImage(
                vector=transform_point_1,
                a=periodic_box.o.a,
                b=periodic_box.o.b,
                c=periodic_box.o.c,
            )
            boolean_math = invert_matrix.o.invertible & (
                periodic_nearest_image.o.vector.length()
                > transform_point_1.length() - 0.0001
            )
            set_position = (
                g.Points(position=transform_point_1, count=17576)
                >> g.DeleteGeometry.point(selection=~boolean_math)
                >> g.MergeByDistance()
                >> g.SetPosition(
                    offset=g.TransformPoint(
                        transform=periodic_box.o.matrix, vector=(0.5, 0.5, 0.5)
                    )
                )
            )
        with g.Frame("Annotations"):
            _string_2 = g.String(
                "The simulation box annotation draws these same cells in Python (draw_triclinic_cell and draw_wigner_seitz_cell in annotations/base.py). Annotations are planned to move to geometry nodes, using this node for the box."
            )
        switch = g.InvertMatrix(periodic_box.o.matrix).o.invertible.switch.geometry(
            true=g.MenuSwitch.geometry(wrap, {"Box": points, "Compact": set_position})
            >> g.ConvexHull()
        )
        with g.Frame("Edges"):
            _string_3 = g.String(
                "Convex Hull triangulates the faces, so keep only the edges between faces that are not coplanar."
            )
            mesh_to_curve = g.MeshToCurve(
                mesh=switch, selection=g.EdgeAngle().o.unsigned_angle > 0.001
            )

        switch >> mesh
        mesh_to_curve >> edges


ASSET = PeriodicCell

ASSET_METADATA = {
    "description": "Mesh and edges of the periodic unit cell, the box spanned by the box vectors or the compact cell around its centre",
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
