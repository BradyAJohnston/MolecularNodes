# Node-group asset 'Selected Instances' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .separate_first_point import SeparateFirstPoint


class SelectedInstances(AssetGeometryGroup):
    """
    Selected Instances

    Parameters
    ----------
    instances : InputGeometry
        Geometry containing instances to check if they are selected or not
    selection : InputBoolean
        The selection to apply to the bounding boxes of the Instances

    Inputs
    ------
    i.instances : GeometrySocket
        Geometry containing instances to check if they are selected or not
    i.selection : BooleanSocket
        The selection to apply to the bounding boxes of the Instances

    Outputs
    -------
    o.entirely_selected : BooleanSocket
        All points of the instances bounding box are selected
    o.partially_selected : BooleanSocket
        Some points of the Instance's bounding box are selected
    o.not_selected : BooleanSocket
        No points of the Instance's bounding box are selected
    """

    _name = "Selected Instances"
    _asset_name = "Selected Instances"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        instances: GeometrySocket
        """Geometry containing instances to check if they are selected or not"""
        selection: BooleanSocket
        """The selection to apply to the bounding boxes of the Instances"""

    class _Outputs(SocketAccessor):
        entirely_selected: BooleanSocket
        """All points of the instances bounding box are selected"""
        partially_selected: BooleanSocket
        """Some points of the Instance's bounding box are selected"""
        not_selected: BooleanSocket
        """No points of the Instance's bounding box are selected"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        instances: InputGeometry = None,
        selection: InputBoolean = False,
    ):
        super().__init__(**{"Instances": instances, "Selection": selection})

    def _build_group(self, tree):
        instances = tree.inputs.geometry(
            "Instances",
            description="Geometry containing instances to check if they are selected or not",
        )
        selection = tree.inputs.boolean(
            "Selection",
            False,
            description="The selection to apply to the bounding boxes of the Instances",
            hide_value=True,
        )
        entirely_selected = tree.outputs.boolean(
            "Entirely Selected",
            description="All points of the instances bounding box are selected",
        )
        partially_selected = tree.outputs.boolean(
            "Partially Selected",
            description="Some points of the Instance's bounding box are selected",
        )
        not_selected = tree.outputs.boolean(
            "Not Selected",
            description="No points of the Instance's bounding box are selected",
        )

        capture = g.CaptureAttribute.instance(
            geometry=g.BoundingBox(geometry=instances)
        )
        index = capture.items.integer("Index", g.Index())
        accumulate_field = g.AccumulateField.point.integer(selection, index.output)
        index_1 = g.Index()
        capture_1 = g.CaptureAttribute.point(
            geometry=capture.o.geometry
            >> g.RealizeInstances(realize_to_point_domain=True)
        )
        all_points_select = capture_1.items.boolean(
            "All Points Select", g.Compare.integer.equal(accumulate_field.o.total, 8)
        )
        some_points_selected = capture_1.items.boolean(
            "Some Points Selected",
            (accumulate_field.o.total < 8) & (accumulate_field.o.total > 0),
        )
        no_points_selected = capture_1.items.boolean(
            "No points selected", g.Compare.integer.equal(accumulate_field.o.total, 0)
        )
        group = SeparateFirstPoint(
            geometry=capture_1.o.geometry, sort=False, group_id=index.output
        )
        sample_index = g.SampleIndex(
            geometry=group,
            value=all_points_select.output,
            index=index_1,
            data_type="BOOLEAN",
        )
        sample_index_1 = g.SampleIndex(
            geometry=group,
            value=some_points_selected.output,
            index=index_1,
            data_type="BOOLEAN",
        )
        sample_index_2 = g.SampleIndex(
            geometry=group,
            value=no_points_selected.output,
            index=index_1,
            data_type="BOOLEAN",
        )

        sample_index >> entirely_selected
        sample_index_1 >> partially_selected
        sample_index_2 >> not_selected


ASSET = SelectedInstances

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
