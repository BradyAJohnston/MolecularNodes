# Node-group asset 'Visualize Relative Atoms' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputVector,
)
from .index_distance import IndexDistance
from .primitive_arrow import PrimitiveArrow


class VisualizeRelativeAtoms(AssetGeometryGroup):
    """
    Visualize Relative Atoms

    Parameters
    ----------
    atoms : InputGeometry
        Atoms
    selection : InputBoolean
        Selection
    scale : InputFloat
        Scale
    position : InputVector
        Position
    target_index : InputInteger
        Index for the selected point to measure to
    target_position : InputVector
        Target Position

    Inputs
    ------
    i.atoms : GeometrySocket
        Atoms
    i.selection : BooleanSocket
        Selection
    i.scale : FloatSocket
        Scale
    i.position : VectorSocket
        Position
    i.target_index : IntegerSocket
        Index for the selected point to measure to
    i.target_position : VectorSocket
        Target Position

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Visualize Relative Atoms"
    _asset_name = "Visualize Relative Atoms"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atoms"""
        selection: BooleanSocket
        """Selection"""
        scale: FloatSocket
        """Scale"""
        position: VectorSocket
        """Position"""
        target_index: IntegerSocket
        """Index for the selected point to measure to"""
        target_position: VectorSocket
        """Target Position"""

    class _Outputs(SocketAccessor):
        instances: GeometrySocket
        """Instances"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        scale: InputFloat = 1.0,
        position: InputVector = None,
        target_index: InputInteger = 100,
        target_position: InputVector = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Scale": scale,
                "Position": position,
                "Target Index": target_index,
                "Target Position": target_position,
            }
        )

    def _build_group(self, tree):
        atoms = tree.inputs.geometry("Atoms")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        scale = tree.inputs.float("Scale", 1.0, min_value=-10_000.0, max_value=10_000.0)
        position = tree.inputs.vector(
            "Position", (0.0, 0.0, 0.0), hide_value=True, default_input="POSITION"
        )
        target_index = tree.inputs.integer(
            "Target Index",
            100,
            description="Index for the selected point to measure to",
            min_value=0,
        )
        target_position = tree.inputs.vector(
            "Target Position",
            (0.0, 0.0, 0.0),
            hide_value=True,
            default_input="POSITION",
        )
        instances = tree.outputs.geometry("Instances")

        group = IndexDistance(target_index=target_index, position=target_position)
        (
            atoms
            >> g.SetPosition(position=position)
            >> g.InstanceOnPoints(
                selection=selection,
                instance=PrimitiveArrow(
                    vertices=3, ratio=0.5625, value=(0.0, 0.0, 0.0, 1.0)
                ),
                rotation=group.o.rotation,
                scale=g.CombineXYZ(z=group.o.distance * scale, x=0.02, y=0.02),
            )
            >> instances
        )


ASSET = VisualizeRelativeAtoms

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
