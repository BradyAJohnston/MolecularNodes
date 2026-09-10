# Node-group asset 'Sample Mixed Rotation' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputGeometry, InputRotation
from .index_mix_rotation import IndexMixRotation


class SampleMixedRotation(AssetGeometryGroup):
    """
    Sample Mixed Rotation

    Parameters
    ----------
    geometry : InputGeometry
        The geometry to sample the values from
    rotation : InputRotation
        The field to mix and evaluate on the sample geometry
    index : InputFloat
        The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes

    Inputs
    ------
    i.geometry : GeometrySocket
        The geometry to sample the values from
    i.rotation : RotationSocket
        The field to mix and evaluate on the sample geometry
    i.index : FloatSocket
        The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes

    Outputs
    -------
    o.rotation : RotationSocket
        The evaluated and mixed field, sampled from the sample geometry at the given `Index`
    """

    _name = "Sample Mixed Rotation"
    _asset_name = "Sample Mixed Rotation"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.sample_mixed_rotation"}

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """The geometry to sample the values from"""
        rotation: RotationSocket
        """The field to mix and evaluate on the sample geometry"""
        index: FloatSocket
        """The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes"""

    class _Outputs(SocketAccessor):
        rotation: RotationSocket
        """The evaluated and mixed field, sampled from the sample geometry at the given `Index`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        rotation: InputRotation = None,
        index: InputFloat = 0.0,
    ):
        super().__init__(**{"Geometry": geometry, "Rotation": rotation, "Index": index})

    def _build_group(self, tree):
        geometry = tree.inputs.geometry(
            "Geometry", description="The geometry to sample the values from"
        )
        rotation = tree.inputs.rotation(
            "Rotation",
            (0.0, 0.0, 0.0),
            description="The field to mix and evaluate on the sample geometry",
            hide_value=True,
        )
        index = tree.inputs.float(
            "Index",
            0.0,
            description="The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes",
        )
        rotation_1 = tree.outputs.rotation(
            "Rotation",
            description="The evaluated and mixed field, sampled from the sample geometry at the given `Index`",
        )

        group = IndexMixRotation(rotation=rotation, index=index)
        (
            geometry
            >> g.SampleIndex(
                value=group.o.rotation,
                index=group.o.from_,
                data_type="QUATERNION",
                clamp=True,
            )
            >> rotation_1
        )


ASSET = SampleMixedRotation

ASSET_METADATA = {
    "catalog_id": "dd5f0199-fa8b-4b01-a972-2dc586a3e60f",
}
