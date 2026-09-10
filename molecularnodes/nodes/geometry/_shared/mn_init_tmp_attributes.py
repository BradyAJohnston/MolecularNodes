# Node group '.MN_init_tmp_attributes' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, GeometrySocket, SocketAccessor
from nodebpy.types import InputGeometry
from ..attribute_run import AttributeRun
from ..group_parameter import GroupParameter


class MN_init_tmp_attributes(CustomGeometryGroup):
    """
    .MN_init_tmp_attributes

    Parameters
    ----------
    geometry : InputGeometry
        Geometry

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = ".MN_init_tmp_attributes"
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry._mn_init_tmp_attributes"}

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

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
    ):
        super().__init__(**{"Geometry": geometry})

    def _build_group(self, tree):
        geometry = tree.inputs.geometry("Geometry")
        geometry_1 = tree.outputs.geometry("Geometry")

        group = GroupParameter(
            group_id=g.NamedAttribute.integer("tmp_ss_ID").o.attribute
        )
        store_named_attribute = g.StoreNamedAttribute.point.integer(
            geometry,
            name="tmp_ss_ID",
            value=AttributeRun(
                name="sec_struct", group_id=g.CurveOfPoint().o.curve_index
            ).o.group_id,
        )
        capture = g.CaptureAttribute.point(geometry=store_named_attribute)
        is_first = capture.items.boolean("Is First", group.o.is_first)
        is_last = capture.items.boolean("Is Last", group.o.is_last)
        size = capture.items.integer("Size", group.o.group_size)
        (
            capture.o.geometry
            >> g.StoreNamedAttribute.point.integer(
                name="tmp_ss_size", value=size.output
            )
            >> g.StoreNamedAttribute.point.boolean(
                name="tmp_ss_first", value=is_first.output
            )
            >> g.StoreNamedAttribute.point.boolean(
                name="tmp_ss_last", value=is_last.output
            )
            >> g.StoreNamedAttribute.point.integer(
                name="tmp_idx_curve", value=g.CurveOfPoint().o.curve_index
            )
            >> geometry_1
        )
