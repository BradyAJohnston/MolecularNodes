# Node-group asset 'Slice Edge Instances' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
)
from nodebpy.types import InputBoolean, InputGeometry
from .selected_instances import SelectedInstances


class SliceEdgeInstances(AssetGeometryGroup):
    """
    Slice Edge Instances

    Parameters
    ----------
    instances : InputGeometry
        Instances
    selection : InputBoolean
        Selection

    Inputs
    ------
    i.instances : GeometrySocket
        Instances
    i.selection : BooleanSocket
        Selection

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    o.realized_points : GeometrySocket
        Realized Points
    """

    _name = "Slice Edge Instances"
    _asset_name = "Slice Edge Instances"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        instances: GeometrySocket
        """Instances"""
        selection: BooleanSocket
        """Selection"""

    class _Outputs(SocketAccessor):
        instances: GeometrySocket
        """Instances"""
        realized_points: GeometrySocket
        """Realized Points"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        instances: InputGeometry = None,
        selection: InputBoolean = True,
    ):
        super().__init__(**{"Instances": instances, "Selection": selection})

    def _build_group(self, tree):
        instances = tree.inputs.geometry("Instances")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        instances_1 = tree.outputs.geometry("Instances")
        realized_points = tree.outputs.geometry("Realized Points")

        group = SelectedInstances(instances=instances, selection=selection)
        separate_geometry = g.SeparateGeometry.instance(
            instances, group.o.entirely_selected
        )
        (
            g.SeparateGeometry.instance(instances, group.o.partially_selected)
            >> g.RealizeInstances(realize_to_point_domain=True)
            >> g.SeparateGeometry.point(selection=selection)
            >> realized_points
        )

        separate_geometry >> instances_1


ASSET = SliceEdgeInstances

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
