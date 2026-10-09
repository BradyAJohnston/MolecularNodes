# Node-group asset "Screen Space Transform" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    MatrixSocket,
    ObjectSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputObject


class ScreenSpaceTransform(AssetGeometryGroup):
    """
    Transform from screen space to the object's local space, on a plane at a set distance in front of the camera. (0, 0) is the bottom left of the camera frame and (1, 1) the top right. One unit is the width and height of the frame, so geometry keeps a constant size on screen and faces the camera

    Parameters
    ----------
    camera : InputObject
        Camera that defines the screen. When empty the scene's active camera is used
    distance : InputFloat
        Distance in front of the camera to place the screen plane, in world units
    normalize : InputBoolean
        Scale X by the frame height instead of the frame width, so X runs from 0 to the aspect ratio and geometry keeps its proportions

    Inputs
    ------
    i.camera : ObjectSocket
        Camera that defines the screen. When empty the scene's active camera is used
    i.distance : FloatSocket
        Distance in front of the camera to place the screen plane, in world units
    i.normalize : BooleanSocket
        Scale X by the frame height instead of the frame width, so X runs from 0 to the aspect ratio and geometry keeps its proportions

    Outputs
    -------
    o.transform : MatrixSocket
        Transform from screen space to the object's local space
    o.aspect : FloatSocket
        Width of the camera frame divided by its height
    """

    _name = "Screen Space Transform"
    _asset_name = "Screen Space Transform"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Transform from screen space to the object's local space, on a plane at a set distance in front of the camera. (0, 0) is the bottom left of the camera frame and (1, 1) the top right. One unit is the width and height of the frame, so geometry keeps a constant size on screen and faces the camera"
    }

    class _Inputs(SocketAccessor):
        camera: ObjectSocket
        """Camera that defines the screen. When empty the scene's active camera is used"""
        distance: FloatSocket
        """Distance in front of the camera to place the screen plane, in world units"""
        normalize: BooleanSocket
        """Scale X by the frame height instead of the frame width, so X runs from 0 to the aspect ratio and geometry keeps its proportions"""

    class _Outputs(SocketAccessor):
        transform: MatrixSocket
        """Transform from screen space to the object's local space"""
        aspect: FloatSocket
        """Width of the camera frame divided by its height"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        camera: InputObject = None,
        distance: InputFloat = 1.0,
        normalize: InputBoolean = False,
    ):
        super().__init__(
            **{"Camera": camera, "Distance": distance, "Normalize": normalize}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        camera = tree.inputs.object(
            "Camera",
            description="Camera that defines the screen. When empty the scene's active camera is used",
            optional_label=True,
        )
        distance = tree.inputs.float(
            "Distance",
            1.0,
            description="Distance in front of the camera to place the screen plane, in world units",
            min_value=0.0,
            subtype="DISTANCE",
        )
        normalize = tree.inputs.boolean(
            "Normalize",
            False,
            description="Scale X by the frame height instead of the frame width, so X runs from 0 to the aspect ratio and geometry keeps its proportions",
        )
        transform = tree.outputs.matrix(
            "Transform",
            description="Transform from screen space to the object's local space",
        )
        aspect = tree.outputs.float(
            "Aspect", description="Width of the camera frame divided by its height"
        )

        with g.Frame("Camera fallback"):
            _string = g.String(
                string="An empty Camera input reports a focal length of 0, in which case the scene's active camera is used instead."
            )
            switch = (g.CameraInfo(camera=camera).o.focal_length > 0.0).switch.object(
                g.ActiveCamera(), camera
            )
        with g.Frame("Unproject frame corners"):
            _string_1 = g.String(
                string="The bottom left and top right corners of the frame are unprojected through the inverse projection matrix at the near and far clip planes, then interpolated along that ray to the plane at Distance in front of the camera. This works for perspective and orthographic cameras."
            )
            invert_matrix = g.CameraInfo(camera=switch).o.projection_matrix.invert()
            project_point = g.ProjectPoint(
                transform=invert_matrix, vector=(-1.0, -1.0, -1.0)
            )
            project_point_1 = g.ProjectPoint(
                transform=invert_matrix, vector=(-1.0, -1.0, 1.0)
            )
            project_point_2 = g.ProjectPoint(
                transform=invert_matrix, vector=(1.0, 1.0, -1.0)
            )
            project_point_3 = g.ProjectPoint(
                transform=invert_matrix, vector=(1.0, 1.0, 1.0)
            )
            math_1 = (distance * -1.0 - project_point.o.vector.z) / (
                project_point_1.o.vector.z - project_point.o.vector.z
            )
            math_2 = (distance * -1.0 - project_point_2.o.vector.z) / (
                project_point_3.o.vector.z - project_point_2.o.vector.z
            )
            vector_math = (
                project_point.o.vector
                + (project_point_1.o.vector - project_point) * math_1
            )
            vector_math_1 = (
                project_point_2.o.vector
                + (project_point_3.o.vector - project_point_2) * math_2
                - vector_math
            )
        with g.Frame("Screen to object space"):
            _string_2 = g.String(
                string="Screen space is scaled to the frame and moved to its corner in camera space, then into world space with the camera's location and rotation (ignoring its scale), then into the local space of the object this node is evaluated on."
            )
            object_info = g.ObjectInfo(object=switch)
            combine_xyz = g.CombineXYZ(
                x=normalize.switch.float(vector_math_1.x, vector_math_1.y),
                y=vector_math_1.y,
                z=vector_math_1.y,
            )
            multiply_matrices = g.MultiplyMatrices(
                matrix=g.CombineTransform(
                    translation=object_info.o.location, rotation=object_info.o.rotation
                ),
                matrix_001=g.CombineTransform(
                    translation=vector_math, scale=combine_xyz
                ),
            )
            multiply_matrices_1 = g.MultiplyMatrices(
                matrix=g.ObjectInfo(object=g.SelfObject()).o.transform.invert(),
                matrix_001=multiply_matrices,
            )
        vector_math_1.x / vector_math_1.y >> aspect

        multiply_matrices_1 >> transform


ASSET = ScreenSpaceTransform

ASSET_METADATA = {
    "description": "Transform from screen space to the object's local space, on a plane at a set distance in front of the camera. (0, 0) is the bottom left of the camera frame and (1, 1) the top right. One unit is the width and height of the frame, so geometry keeps a constant size on screen and faces the camera",
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
