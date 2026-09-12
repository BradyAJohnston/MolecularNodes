# Node-group asset "Rotation from ZYZ" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
)
from nodebpy.types import InputFloat


class RotationFromZYZ(AssetGeometryGroup):
    """
    Combine a rotation defined as ZYZ common in electron tomography

    Parameters
    ----------
    phi : InputFloat
        First rotation around the Z axis
    theta : InputFloat
        Second rotation around the Y axis
    psi : InputFloat
        Third rotation around the Z axis

    Inputs
    ------
    i.phi : FloatSocket
        First rotation around the Z axis
    i.theta : FloatSocket
        Second rotation around the Y axis
    i.psi : FloatSocket
        Third rotation around the Z axis

    Outputs
    -------
    o.rotation : RotationSocket
        The combined Rotation
    """

    _name = "Rotation from ZYZ"
    _asset_name = "Rotation from ZYZ"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Combine a rotation defined as ZYZ common in electron tomography",
        "node_tool_idname": "geometry.rotation_from_zyz",
    }

    class _Inputs(SocketAccessor):
        phi: FloatSocket
        """First rotation around the Z axis"""
        theta: FloatSocket
        """Second rotation around the Y axis"""
        psi: FloatSocket
        """Third rotation around the Z axis"""

    class _Outputs(SocketAccessor):
        rotation: RotationSocket
        """The combined Rotation"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        phi: InputFloat = 0.0,
        theta: InputFloat = 0.0,
        psi: InputFloat = 0.0,
    ):
        super().__init__(**{"Phi": phi, "Theta": theta, "Psi": psi})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        phi = tree.inputs.float(
            "Phi",
            0.0,
            description="First rotation around the Z axis",
            min_value=-10_000.0,
            max_value=10_000.0,
            subtype="ANGLE",
        )
        theta = tree.inputs.float(
            "Theta",
            0.0,
            description="Second rotation around the Y axis",
            min_value=-10_000.0,
            max_value=10_000.0,
            subtype="ANGLE",
        )
        psi = tree.inputs.float(
            "Psi",
            0.0,
            description="Third rotation around the Z axis",
            min_value=-10_000.0,
            max_value=10_000.0,
            subtype="ANGLE",
        )
        rotation = tree.outputs.rotation(
            "Rotation", description="The combined Rotation"
        )

        rotate_rotation = g.AxisAngleToRotation(angle=phi * -1.0).o.rotation.rotate(
            g.AxisAngleToRotation(angle=theta * -1.0, axis=(0.0, 1.0, 0.0))
        )
        rotate_rotation.rotate(g.AxisAngleToRotation(angle=psi * -1.0)) >> rotation


ASSET = RotationFromZYZ

ASSET_METADATA = {
    "description": "Combine a rotation defined as ZYZ common in electron tomography",
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
