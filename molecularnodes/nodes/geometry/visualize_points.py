# Node-group asset 'Visualize Points' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputGeometry, InputVector
from .primitive_arrow import PrimitiveArrow
from .vector_from_point import VectorFromPoint


class VisualizePoints(AssetGeometryGroup):
    """
    Visualize Points

    Parameters
    ----------
    points : InputGeometry
        Points
    selection : InputBoolean
        Selection
    target : InputVector
        Vector that is the target
    position : InputVector
        Position of the current point

    Inputs
    ------
    i.points : GeometrySocket
        Points
    i.selection : BooleanSocket
        Selection
    i.target : VectorSocket
        Vector that is the target
    i.position : VectorSocket
        Position of the current point

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Visualize Points"
    _asset_name = "Visualize Points"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        points: GeometrySocket
        """Points"""
        selection: BooleanSocket
        """Selection"""
        target: VectorSocket
        """Vector that is the target"""
        position: VectorSocket
        """Position of the current point"""

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
        points: InputGeometry = None,
        selection: InputBoolean = True,
        target: InputVector = None,
        position: InputVector = None,
    ):
        super().__init__(
            **{
                "Points": points,
                "Selection": selection,
                "Target": target,
                "Position": position,
            }
        )

    def _build_group(self, tree):
        points = tree.inputs.geometry("Points")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        target = tree.inputs.vector(
            "Target",
            (0.0, 0.0, 0.0),
            description="Vector that is the target",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        position = tree.inputs.vector(
            "Position",
            (0.0, 0.0, 0.0),
            description="Position of the current point",
            min_value=-10_000.0,
            max_value=10_000.0,
            hide_value=True,
            default_input="POSITION",
        )
        instances = tree.outputs.geometry("Instances")

        group = VectorFromPoint(target=target, position=position)
        (
            points
            >> g.InstanceOnPoints(
                selection=selection,
                instance=PrimitiveArrow(value=(0.0, 0.0, 0.0, 1.0)),
                rotation=group.o.rotation,
                scale=g.CombineXYZ(z=group.o.length, x=0.5, y=0.5),
            )
            >> instances
        )


ASSET = VisualizePoints

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
