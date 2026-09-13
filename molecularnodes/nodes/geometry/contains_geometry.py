# Node-group asset "Contains Geometry" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputGeometry


class ContainsGeometry(AssetGeometryGroup):
    """
    Contains Geometry

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
    o.not_empty : BooleanSocket
        True when the input contains some geometry
    o.empty : BooleanSocket
        True when the input contains no geometry
    """

    _name = "Contains Geometry"
    _asset_name = "Contains Geometry"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

    class _Outputs(SocketAccessor):
        not_empty: BooleanSocket
        """True when the input contains some geometry"""
        empty: BooleanSocket
        """True when the input contains no geometry"""

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

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        not_empty = tree.outputs.boolean(
            "Not Empty", description="True when the input contains some geometry"
        )
        empty = tree.outputs.boolean(
            "Empty", description="True when the input contains no geometry"
        )

        compare = g.Compare.integer.equal(
            g.DomainSize(geometry=geometry, component="POINTCLOUD").o.point_count, 0
        )
        boolean_math = (
            g.Compare.integer.equal(
                g.DomainSize(geometry=geometry).o.point_count, 0
            ).o.result
            & compare
        )
        compare_1 = g.Compare.integer.equal(
            g.DomainSize(geometry=geometry, component="CURVE").o.point_count, 0
        )
        compare_2 = g.Compare.integer.equal(
            g.DomainSize(geometry=geometry, component="INSTANCES").o.instance_count, 0
        )
        compare_3 = g.Compare.integer.equal(
            g.DomainSize(geometry=geometry, component="GREASEPENCIL").o.layer_count, 0
        )
        boolean_math_1 = boolean_math & compare_1 & compare_2 & compare_3
        ~boolean_math_1 >> not_empty

        boolean_math_1 >> empty


ASSET = ContainsGeometry

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
