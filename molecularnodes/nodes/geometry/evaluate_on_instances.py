# Node-group asset 'Evaluate on Instances' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ClosureSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputClosure, InputGeometry
from .evaluate_per_group import EvaluatePerGroup


class EvaluateOnInstances(AssetGeometryGroup):
    """
    Evaluate on Instances

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to split into two parts
    closure : InputClosure
        Closure

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to split into two parts
    i.closure : ClosureSocket
        Closure

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Evaluate on Instances"
    _asset_name = "Evaluate on Instances"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to split into two parts"""
        closure: ClosureSocket
        """Closure"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        closure: InputClosure = None,
    ):
        super().__init__(**{"Geometry": geometry, "Closure": closure})

    def _build_group(self, tree):
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to split into two parts"
        )
        closure = tree.inputs.closure("Closure")
        geometry_1 = tree.outputs.geometry("Geometry")

        with g.Frame("Get unique geometry references and evaluate on them"):
            instance_reference = g.InstanceReference()
            separate_geometry = g.SeparateGeometry.instance(
                geometry,
                ~g.AccumulateField.instance.integer(
                    group_index=instance_reference
                ).o.trailing,
            )
            capture = g.CaptureAttribute.instance(
                geometry=g.SortElements.instance(
                    separate_geometry.o.selection, sort_weight=instance_reference
                )
            )
            index = capture.items.integer("Index", g.Index())
            with g.Frame("Clear Instance Transforms"):
                realize_instances = (
                    capture.o.geometry
                    >> g.SetInstanceTransform(transform=g.CombineMatrix())
                    >> g.RealizeInstances()
                )
            group = EvaluatePerGroup(
                geometry=realize_instances,
                closure=closure,
                group="Group ID",
                group_id=index.output,
            )
        index_1 = g.Index()
        sample_index = g.SampleIndex(
            geometry=geometry,
            value=g.InstanceTransform(),
            index=index_1,
            data_type="FLOAT4X4",
            domain="INSTANCE",
        )
        sample_index_1 = g.SampleIndex(
            geometry=geometry,
            value=g.InstanceReference(),
            index=index_1,
            data_type="INT",
            domain="INSTANCE",
        )
        (
            g.InstancesToPoints(instances=geometry, radius=0.05)
            >> g.InstanceOnPoints(
                instance=group.o.instances,
                instance_index=sample_index_1,
                pick_instance=True,
            )
            >> g.SetInstanceTransform(transform=sample_index)
            >> geometry_1
        )


ASSET = EvaluateOnInstances

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
