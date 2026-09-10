# Node group 'Smooth by Angle' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry


class SmoothByAngle(CustomGeometryGroup):
    """
    Set the sharpness of mesh edges based on the angle between the neighboring faces

    Parameters
    ----------
    mesh : InputGeometry
        Mesh
    angle : InputFloat
        Maximum face angle for smooth edges
    ignore_sharpness : InputBoolean
        Ignore Sharpness

    Inputs
    ------
    i.mesh : GeometrySocket
        Mesh
    i.angle : FloatSocket
        Maximum face angle for smooth edges
    i.ignore_sharpness : BooleanSocket
        Ignore Sharpness

    Outputs
    -------
    o.mesh : GeometrySocket
        Mesh
    """

    _name = "Smooth by Angle"
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Set the sharpness of mesh edges based on the angle between the neighboring faces",
        "node_tool_idname": "geometry.smooth_by_angle",
        "is_modifier": True,
    }

    class _Inputs(SocketAccessor):
        mesh: GeometrySocket
        """Mesh"""
        angle: FloatSocket
        """Maximum face angle for smooth edges"""
        ignore_sharpness: BooleanSocket
        """Ignore Sharpness"""

    class _Outputs(SocketAccessor):
        mesh: GeometrySocket
        """Mesh"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        mesh: InputGeometry = None,
        angle: InputFloat = math.pi / 6,
        ignore_sharpness: InputBoolean = False,
    ):
        super().__init__(
            **{"Mesh": mesh, "Angle": angle, "Ignore Sharpness": ignore_sharpness}
        )

    def _build_group(self, tree):
        mesh = tree.inputs.geometry("Mesh")
        angle = tree.inputs.float(
            "Angle",
            math.pi / 6,
            description="Maximum face angle for smooth edges",
            min_value=0.0,
            max_value=math.pi,
            subtype="ANGLE",
        )
        ignore_sharpness = tree.inputs.boolean(
            "Ignore Sharpness", False, structure_type="SINGLE", force_non_field=True
        )
        mesh_1 = tree.outputs.geometry("Mesh")

        boolean_math = (g.EdgeAngle().o.unsigned_angle <= angle) & (
            g.IsFaceSmooth().o.smooth | ignore_sharpness
        )
        (
            mesh
            >> g.SetShadeSmooth.edge(
                selection=g.IsEdgeSmooth().o.smooth | ignore_sharpness,
                shade_smooth=boolean_math,
            )
            >> g.SetShadeSmooth.face()
            >> mesh_1
        )
