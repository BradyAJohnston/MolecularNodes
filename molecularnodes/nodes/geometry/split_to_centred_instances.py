# Node-group asset "Split to Centred Instances" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputGeometry, InputInteger
from .centroid import Centroid
from .separate_first_point import SeparateFirstPoint


class SplitToCentredInstances(AssetGeometryGroup):
    """
    Split to Centred Instances

    Parameters
    ----------
    geometry : InputGeometry
        The input geometry containing points
    selection : InputBoolean
        Selected points will be used to calculate centre point for Instance
    group_id : InputInteger
        The `Group ID` which will determine how the points are split into their different instances

    Inputs
    ------
    i.geometry : GeometrySocket
        The input geometry containing points
    i.selection : BooleanSocket
        Selected points will be used to calculate centre point for Instance
    i.group_id : IntegerSocket
        The `Group ID` which will determine how the points are split into their different instances

    Outputs
    -------
    o.instances : GeometrySocket
        The points that have been split into their instances
    """

    _name = "Split to Centred Instances"
    _asset_name = "Split to Centred Instances"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.split_to_centred_instances"}

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """The input geometry containing points"""
        selection: BooleanSocket
        """Selected points will be used to calculate centre point for Instance"""
        group_id: IntegerSocket
        """The `Group ID` which will determine how the points are split into their different instances"""

    class _Outputs(SocketAccessor):
        instances: GeometrySocket
        """The points that have been split into their instances"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        selection: InputBoolean = True,
        group_id: InputInteger = 0,
    ):
        super().__init__(
            **{"Geometry": geometry, "Selection": selection, "Group ID": group_id}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry",
            description="The input geometry containing points",
            hide_value=True,
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selected points will be used to calculate centre point for Instance",
            hide_value=True,
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="The `Group ID` which will determine how the points are split into their different instances",
            hide_value=True,
        )
        instances = tree.outputs.geometry(
            "Instances",
            description="The points that have been split into their instances",
        )

        capture = g.CaptureAttribute.point(geometry=geometry)
        centroid = capture.items.vector(
            "Centroid", Centroid(selection=selection, group_id=group_id)
        )
        with g.Frame("Get first point in each group, to be sampled by instance"):
            sample_index = SeparateFirstPoint(
                geometry=capture.o.geometry, group_id=group_id
            ) >> g.SampleIndex(
                value=centroid.output, index=g.Index(), data_type="FLOAT_VECTOR"
            )
        with g.Frame("Centres each group on world origin"):
            set_position = g.SetPosition(
                geometry=capture.o.geometry, offset=centroid.output * -1.0
            )
        with g.Frame("Split to Instances and return to original positon"):
            (
                set_position
                >> g.SplitToInstances.point(group_id=group_id)
                >> g.SetPosition(offset=sample_index)
                >> instances
            )


ASSET = SplitToCentredInstances

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
