# Node group "Project with Depth" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    MatrixSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputFloat, InputMatrix, InputVector


class ProjectWithDepth(CustomGeometryGroup):
    """
    Unproject normalized screen coordinates at a depth from the camera through the inverse projection matrix, for perspective and orthographic cameras. Adapted from the Project with Depth group in Blender's Essentials asset library, which used the perspective depth mapping for orthographic cameras too

    Parameters
    ----------
    normalized : InputVector
        Normalized
    depth : InputFloat
        Depth
    projection : InputMatrix
        Projection
    transform : InputMatrix
        Transform
    clip_start : InputFloat
        Clip Start
    clip_end : InputFloat
        Clip End
    is_orthographic : InputBoolean
        Map Depth to normalized device coordinates linearly, as an orthographic camera does, instead of with the perspective formula

    Inputs
    ------
    i.normalized : VectorSocket
        Normalized
    i.depth : FloatSocket
        Depth
    i.projection : MatrixSocket
        Projection
    i.transform : MatrixSocket
        Transform
    i.clip_start : FloatSocket
        Clip Start
    i.clip_end : FloatSocket
        Clip End
    i.is_orthographic : BooleanSocket
        Map Depth to normalized device coordinates linearly, as an orthographic camera does, instead of with the perspective formula

    Outputs
    -------
    o.vector : VectorSocket
        Vector
    """

    _name = "Project with Depth"
    _color_tag = "VECTOR"
    _tree_properties = {
        "description": "Unproject normalized screen coordinates at a depth from the camera through the inverse projection matrix, for perspective and orthographic cameras. Adapted from the Project with Depth group in Blender's Essentials asset library, which used the perspective depth mapping for orthographic cameras too"
    }

    class _Inputs(SocketAccessor):
        normalized: VectorSocket
        """Normalized"""
        depth: FloatSocket
        """Depth"""
        projection: MatrixSocket
        """Projection"""
        transform: MatrixSocket
        """Transform"""
        clip_start: FloatSocket
        """Clip Start"""
        clip_end: FloatSocket
        """Clip End"""
        is_orthographic: BooleanSocket
        """Map Depth to normalized device coordinates linearly, as an orthographic camera does, instead of with the perspective formula"""

    class _Outputs(SocketAccessor):
        vector: VectorSocket
        """Vector"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        normalized: InputVector = None,
        depth: InputFloat = 1.0,
        projection: InputMatrix = None,
        transform: InputMatrix = None,
        clip_start: InputFloat = 0.1,
        clip_end: InputFloat = 100.0,
        is_orthographic: InputBoolean = False,
    ):
        super().__init__(
            **{
                "Normalized": normalized,
                "Depth": depth,
                "Projection": projection,
                "Transform": transform,
                "Clip Start": clip_start,
                "Clip End": clip_end,
                "Is Orthographic": is_orthographic,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        normalized = tree.inputs.vector(
            "Normalized", (0.5, 0.5), dimensions=2, min_value=0.0, max_value=1.0
        )
        depth = tree.inputs.float(
            "Depth", 1.0, min_value=-10_000.0, max_value=10_000.0, subtype="DISTANCE"
        )
        projection = tree.inputs.matrix("Projection")
        transform = tree.inputs.matrix("Transform")
        clip_start = tree.inputs.float(
            "Clip Start",
            0.1,
            min_value=-10_000.0,
            max_value=10_000.0,
            subtype="DISTANCE",
        )
        clip_end = tree.inputs.float(
            "Clip End",
            100.0,
            min_value=-10_000.0,
            max_value=10_000.0,
            subtype="DISTANCE",
        )
        is_orthographic = tree.inputs.boolean(
            "Is Orthographic",
            False,
            description="Map Depth to normalized device coordinates linearly, as an orthographic camera does, instead of with the perspective formula",
        )
        vector = tree.outputs.vector("Vector", subtype="XYZ")

        with g.Frame("Project Depth"):
            with g.Frame("(f + n) / (f - n)"):
                math_1 = (clip_start + clip_end) / (clip_end - clip_start)
            with g.Frame("2fn / Z*(f - n)"):
                math_2 = clip_start * clip_end * 2.0 / (depth * (clip_end - clip_start))
            math_3 = math_1 - math_2
        with g.Frame("Orthographic depth"):
            _string = g.String(
                string="An orthographic projection maps depth linearly: (2Z - f - n) / (f - n)."
            )
            switch = is_orthographic.switch.float(
                math_3, (depth * 2.0 - clip_start - clip_end) / (clip_end - clip_start)
            )
        with g.Frame():
            vector_1 = normalized.mul_add((2.0, 2.0, 2.0), (-1.0, -1.0, -1.0))
            (
                g.CombineXYZ(x=vector_1.x, y=vector_1.y, z=switch)
                .o.vector.project_point(projection)
                .transform(transform)
                >> vector
            )
