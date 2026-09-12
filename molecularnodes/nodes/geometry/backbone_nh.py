# Node-group asset "Backbone NH" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputMenu
from .angstrom_to_world import AngstromToWorld
from .backbone_c import BackboneC
from .backbone_ca import BackboneCA
from .backbone_n import BackboneN
from .fallback_vector import FallbackVector
from .offset_vector import OffsetVector
from .vector_direction import VectorDirection


class BackboneNH(AssetGeometryGroup):
    """
    Backbone NH

    Parameters
    ----------
    menu : InputMenu | Literal["Read", "Compute"]
        Menu

    Inputs
    ------
    i.menu : MenuSocket
        Menu

    Outputs
    -------
    o.nh : VectorSocket
        NH
    """

    _name = "Backbone NH"
    _asset_name = "Backbone NH"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.backbone_nh"}

    class _Inputs(SocketAccessor):
        menu: MenuSocket
        """Menu"""

    class _Outputs(SocketAccessor):
        nh: VectorSocket
        """NH"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        menu: InputMenu | Literal["Read", "Compute"] = "Compute",
    ):
        super().__init__(**{"Menu": menu})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        menu = tree.inputs.menu("Menu", expanded=True, optional_label=True)
        nh = tree.outputs.vector("NH")

        group = BackboneN(method="Read")
        group_1 = BackboneN()
        mix = g.Mix(
            a_vector=BackboneCA(),
            b_vector=OffsetVector(vector=BackboneC(), offset=-1),
            factor_float=0.5,
            data_type="VECTOR",
            clamp_factor=True,
        )
        string = g.String(string="backbone_NH")
        vector_math = (
            VectorDirection(to=group, from_=BackboneCA(method="Read")).o.direction
            + VectorDirection(to=group, from_=BackboneC(method="Read")).o.direction
        )
        _vector_math_1 = g.VectorMath.multiply_add(
            vector_math.normalize(),
            AngstromToWorld(angstrom=1.01),
            BackboneN(method="Read"),
        )
        vector_math_2 = g.VectorMath.multiply_add(
            VectorDirection(to=group_1, from_=mix.o.result_vector).o.direction,
            AngstromToWorld(angstrom=1.01),
            group_1,
        )
        (
            g.MenuSwitch.vector(
                menu,
                {
                    "Read": g.NamedAttribute.vector(string).o.attribute,
                    "Compute": FallbackVector(
                        name=string, fallback=vector_math_2.o.vector
                    ),
                },
            )
            >> nh
        )

        menu.default_value = "Compute"


ASSET = BackboneNH

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
