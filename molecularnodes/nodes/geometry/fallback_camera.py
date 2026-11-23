# Node-group asset "Fallback Camera" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    MatrixSocket,
    ObjectSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputObject


class FallbackCamera(AssetGeometryGroup):
    """
    Fallback Camera

    Parameters
    ----------
    camera : InputObject
        Camera

    Inputs
    ------
    i.camera : ObjectSocket
        Camera

    Outputs
    -------
    o.transform : MatrixSocket
        Transformation matrix containing the location, rotation and scale of the object
    o.camera : ObjectSocket
        Camera
    """

    _name = "Fallback Camera"
    _asset_name = "Fallback Camera"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        camera: ObjectSocket
        """Camera"""

    class _Outputs(SocketAccessor):
        transform: MatrixSocket
        """Transformation matrix containing the location, rotation and scale of the object"""
        camera: ObjectSocket
        """Camera"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        camera: InputObject = None,
    ):
        super().__init__(**{"Camera": camera})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        camera = tree.inputs.object("Camera")
        transform = tree.outputs.matrix(
            "Transform",
            description="Transformation matrix containing the location, rotation and scale of the object",
        )
        camera_1 = tree.outputs.object("Camera")

        with g.Frame("Fallback for the active camera"):
            object_info = g.ObjectInfo(object=camera, transform_space="RELATIVE")
            compare = g.Compare(
                a=object_info.o.location,
                b=(0.0, 0.0, 0.0),
                epsilon=0.001,
                operation="EQUAL",
                data_type="VECTOR",
                mode="ELEMENT",
            )
            compare_1 = g.Compare(
                a=object_info.o.scale,
                b=(0.0, 0.0, 0.0),
                epsilon=0.001,
                operation="EQUAL",
                data_type="VECTOR",
                mode="ELEMENT",
            )
            compare_2 = g.Compare(
                a=object_info.o.rotation,
                b=(0.0, 0.0, 0.0),
                epsilon=0.001,
                operation="EQUAL",
                data_type="VECTOR",
                mode="ELEMENT",
            )
            _string = g.String(
                'Until we can select "Active Camera" as the default object input, we instead detect the active camera by seeing if any of the transforms are non-default.'
            )
            switch = (compare.o.result & compare_2 & compare_1).switch.object(
                camera, g.ActiveCamera()
            )
        object_info_1 = g.ObjectInfo(object=switch, transform_space="RELATIVE")

        object_info_1 >> transform
        switch >> camera_1


ASSET = FallbackCamera

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
