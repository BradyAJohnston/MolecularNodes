# Node-group asset "Screen Space Geometry" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    GeometrySocket,
    ObjectSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry, InputObject
from .screen_space_transform import ScreenSpaceTransform


class ScreenSpaceGeometry(AssetGeometryGroup):
    """
    Place geometry laid out in screen space in front of the camera. X and Y from 0 to 1 span the camera frame from bottom left to top right, and the geometry faces the camera at a constant size on screen. Use it for labels and annotations that stay fixed in the frame

    Parameters
    ----------
    geometry : InputGeometry
        Geometry laid out in screen space, where X and Y from 0 to 1 span the camera frame
    camera : InputObject
        Camera that defines the screen. When empty the scene's active camera is used
    distance : InputFloat
        Distance in front of the camera to place the geometry, in world units
    normalize : InputBoolean
        Scale X by the frame height instead of the frame width, so X runs from 0 to the aspect ratio and geometry keeps its proportions

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry laid out in screen space, where X and Y from 0 to 1 span the camera frame
    i.camera : ObjectSocket
        Camera that defines the screen. When empty the scene's active camera is used
    i.distance : FloatSocket
        Distance in front of the camera to place the geometry, in world units
    i.normalize : BooleanSocket
        Scale X by the frame height instead of the frame width, so X runs from 0 to the aspect ratio and geometry keeps its proportions

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry placed in front of the camera
    o.aspect : FloatSocket
        Width of the camera frame divided by its height
    """

    _name = "Screen Space Geometry"
    _asset_name = "Screen Space Geometry"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Place geometry laid out in screen space in front of the camera. X and Y from 0 to 1 span the camera frame from bottom left to top right, and the geometry faces the camera at a constant size on screen. Use it for labels and annotations that stay fixed in the frame"
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry laid out in screen space, where X and Y from 0 to 1 span the camera frame"""
        camera: ObjectSocket
        """Camera that defines the screen. When empty the scene's active camera is used"""
        distance: FloatSocket
        """Distance in front of the camera to place the geometry, in world units"""
        normalize: BooleanSocket
        """Scale X by the frame height instead of the frame width, so X runs from 0 to the aspect ratio and geometry keeps its proportions"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry placed in front of the camera"""
        aspect: FloatSocket
        """Width of the camera frame divided by its height"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        camera: InputObject = None,
        distance: InputFloat = 1.0,
        normalize: InputBoolean = False,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Camera": camera,
                "Distance": distance,
                "Normalize": normalize,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry",
            description="Geometry laid out in screen space, where X and Y from 0 to 1 span the camera frame",
        )
        camera = tree.inputs.object(
            "Camera",
            description="Camera that defines the screen. When empty the scene's active camera is used",
            optional_label=True,
        )
        distance = tree.inputs.float(
            "Distance",
            1.0,
            description="Distance in front of the camera to place the geometry, in world units",
            min_value=0.0,
            subtype="DISTANCE",
        )
        normalize = tree.inputs.boolean(
            "Normalize",
            False,
            description="Scale X by the frame height instead of the frame width, so X runs from 0 to the aspect ratio and geometry keeps its proportions",
        )
        geometry_1 = tree.outputs.geometry(
            "Geometry", description="Geometry placed in front of the camera"
        )
        aspect = tree.outputs.float(
            "Aspect", description="Width of the camera frame divided by its height"
        )

        screen_space_transform = ScreenSpaceTransform(
            camera=camera, distance=distance, normalize=normalize
        )
        (
            geometry
            >> g.TransformGeometry(
                transform=screen_space_transform.o.transform, mode="Matrix"
            )
            >> geometry_1
        )

        screen_space_transform.o.aspect >> aspect


ASSET = ScreenSpaceGeometry

ASSET_METADATA = {
    "description": "Place geometry laid out in screen space in front of the camera. X and Y from 0 to 1 span the camera frame from bottom left to top right, and the geometry faces the camera at a constant size on screen. Use it for labels and annotations that stay fixed in the frame",
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
