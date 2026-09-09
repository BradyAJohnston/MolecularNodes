# Node-group asset 'TEM Rotation' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputString
from .rotation_from_zyz import RotationFromZYZ


class TEMRotation(AssetGeometryGroup):
    """
    TEM Rotation

    Parameters
    ----------
    phi : InputString
        Phi
    theta : InputString
        Theta
    psi : InputString
        Psi

    Inputs
    ------
    i.phi : StringSocket
        Phi
    i.theta : StringSocket
        Theta
    i.psi : StringSocket
        Psi

    Outputs
    -------
    o.rotation : RotationSocket
        The combined Rotation
    o.boolean : BooleanSocket
        Boolean
    """

    _name = "TEM Rotation"
    _asset_name = "TEM Rotation"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        phi: StringSocket
        """Phi"""
        theta: StringSocket
        """Theta"""
        psi: StringSocket
        """Psi"""

    class _Outputs(SocketAccessor):
        rotation: RotationSocket
        """The combined Rotation"""
        boolean: BooleanSocket
        """Boolean"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        phi: InputString = "rlnAngleRot",
        theta: InputString = "rlnAngleTilt",
        psi: InputString = "rlnAnglePsi",
    ):
        super().__init__(**{"Phi": phi, "Theta": theta, "Psi": psi})

    def _build_group(self, tree):
        phi = tree.inputs.string("Phi", "rlnAngleRot", optional_label=True)
        theta = tree.inputs.string("Theta", "rlnAngleTilt", optional_label=True)
        psi = tree.inputs.string("Psi", "rlnAnglePsi", optional_label=True)
        rotation = tree.outputs.rotation(
            "Rotation", description="The combined Rotation"
        )
        boolean = tree.outputs.boolean("Boolean")

        named_attribute = g.NamedAttribute.float(phi)
        named_attribute_1 = g.NamedAttribute.float(theta)
        named_attribute_2 = g.NamedAttribute.float(psi)
        (
            (
                named_attribute.o.exists
                & (named_attribute_1.o.exists & named_attribute_2.o.exists)
            )
            >> boolean
        )
        (
            RotationFromZYZ(
                phi=named_attribute.o.attribute.to_radians(),
                theta=named_attribute_1.o.attribute.to_radians(),
                psi=named_attribute_2.o.attribute.to_radians(),
            )
            >> rotation
        )


ASSET = TEMRotation

ASSET_METADATA = {
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
