# Node group "Screen to 3D Space" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    FloatSocket,
    ObjectSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputObject, InputVector
from .project_with_depth import ProjectWithDepth


class ScreenTo3DSpace(CustomGeometryGroup):
    """
    Project normalized screen coordinates to 3D at a depth from the camera, in the local space of this object. Adapted from the Screen to 3D Space group in Blender's Essentials asset library

    Parameters
    ----------
    normalized : InputVector
        2D screen space coordinates in the [0, 1] range
    depth : InputFloat
        Depth from camera in scene units. Same as depth pass
    camera : InputObject
        The camera used for rendering the scene

    Inputs
    ------
    i.normalized : VectorSocket
        2D screen space coordinates in the [0, 1] range
    i.depth : FloatSocket
        Depth from camera in scene units. Same as depth pass
    i.camera : ObjectSocket
        The camera used for rendering the scene

    Outputs
    -------
    o.vector : VectorSocket
        Projected 3D vector in 3D space
    """

    _name = "Screen to 3D Space"
    _color_tag = "VECTOR"
    _tree_properties = {
        "description": "Project normalized screen coordinates to 3D at a depth from the camera, in the local space of this object. Adapted from the Screen to 3D Space group in Blender's Essentials asset library"
    }

    class _Inputs(SocketAccessor):
        normalized: VectorSocket
        """2D screen space coordinates in the [0, 1] range"""
        depth: FloatSocket
        """Depth from camera in scene units. Same as depth pass"""
        camera: ObjectSocket
        """The camera used for rendering the scene"""

    class _Outputs(SocketAccessor):
        vector: VectorSocket
        """Projected 3D vector in 3D space"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        normalized: InputVector = None,
        depth: InputFloat = 0.5,
        camera: InputObject = None,
    ):
        super().__init__(**{"Normalized": normalized, "Depth": depth, "Camera": camera})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        normalized = tree.inputs.vector(
            "Normalized",
            (0.5, 0.5),
            description="2D screen space coordinates in the [0, 1] range",
            dimensions=2,
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        depth = tree.inputs.float(
            "Depth",
            0.5,
            description="Depth from camera in scene units. Same as depth pass",
            min_value=-10_000.0,
            max_value=10_000.0,
            subtype="DISTANCE",
        )
        camera = tree.inputs.object(
            "Camera",
            description="The camera used for rendering the scene",
            optional_label=True,
        )
        vector = tree.outputs.vector(
            "Vector", description="Projected 3D vector in 3D space", subtype="XYZ"
        )

        camera_info = g.CameraInfo(camera=camera)
        (
            ProjectWithDepth(
                normalized=normalized,
                depth=depth,
                projection=camera_info.o.projection_matrix.invert(),
                transform=g.ObjectInfo(
                    object=camera, transform_space="RELATIVE"
                ).o.transform,
                clip_start=camera_info.o.clip_start,
                clip_end=camera_info.o.clip_end,
                is_orthographic=camera_info.o.is_orthographic,
            )
            >> vector
        )
