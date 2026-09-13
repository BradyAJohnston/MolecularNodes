# Node group "Animate Collection Pick" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CollectionSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputCollection, InputFloat


class AnimateCollectionPick(CustomGeometryGroup):
    """
    Pick items from a collection based on the index given. The current and next items in the collection are given for interpolation

    Parameters
    ----------
    collection : InputCollection
        Collection
    realize_instances : InputBoolean
        Realize Instances
    item : InputFloat
        Item

    Inputs
    ------
    i.collection : CollectionSocket
        Collection
    i.realize_instances : BooleanSocket
        Realize Instances
    i.item : FloatSocket
        Item

    Outputs
    -------
    o.current : GeometrySocket
        Current
    o.next : GeometrySocket
        Next
    """

    _name = "Animate Collection Pick"
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Pick items from a collection based on the index given. The current and next items in the collection are given for interpolation",
        "node_tool_idname": "geometry.animate_collection_pick",
    }

    class _Inputs(SocketAccessor):
        collection: CollectionSocket
        """Collection"""
        realize_instances: BooleanSocket
        """Realize Instances"""
        item: FloatSocket
        """Item"""

    class _Outputs(SocketAccessor):
        current: GeometrySocket
        """Current"""
        next: GeometrySocket
        """Next"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        collection: InputCollection = None,
        realize_instances: InputBoolean = True,
        item: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Collection": collection,
                "Realize Instances": realize_instances,
                "Item": item,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        collection = tree.inputs.collection("Collection", optional_label=True)
        realize_instances = tree.inputs.boolean("Realize Instances", True)
        item = tree.inputs.float("Item", 1.0, min_value=0.0, max_value=10_000.0)
        current = tree.outputs.geometry("Current")
        next = tree.outputs.geometry("Next")

        collection_info = g.CollectionInfo(
            collection=collection, separate_children=True, transform_space="RELATIVE"
        )
        float_to_integer = g.FloatToInteger(float=item, rounding_mode="FLOOR")
        integer_math = (
            g.DomainSize(
                geometry=collection_info, component="INSTANCES"
            ).o.instance_count
            - 1
        )
        compare = g.Compare.integer.equal(
            g.Index(), g.IntegerMath.minimum(float_to_integer, integer_math)
        )
        compare_1 = g.Compare.integer.equal(
            g.IntegerMath.minimum(float_to_integer.o.integer + 1, integer_math),
            g.Index(),
        )
        (
            g.SeparateGeometry.instance(collection_info, compare)
            >> g.RealizeInstances(
                realize_all=realize_instances, realize_to_point_domain=True
            )
            >> current
        )
        (
            g.SeparateGeometry.instance(collection_info, compare_1)
            >> g.RealizeInstances(
                realize_all=realize_instances, realize_to_point_domain=True
            )
            >> next
        )
