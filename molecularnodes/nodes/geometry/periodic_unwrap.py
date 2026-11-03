# Node-group asset "Periodic Unwrap" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    Default,
    InputBoolean,
    InputFloat,
    InputInteger,
    InputMenu,
    InputVector,
)
from ._shared.periodic_nearest_image import PeriodicNearestImage
from .periodic_box import PeriodicBox
from .periodic_offset import PeriodicOffset


class PeriodicUnwrap(AssetGeometryGroup):
    """
    Make groups of atoms whole across the periodic box, such as molecules split by its boundaries

    Parameters
    ----------
    position : InputVector
        Positions to unwrap. When unconnected: The position from the context.
    group_id : InputInteger
        Atoms with the same ID are kept whole, such as `fragment_id` for bonded molecules or `chain_id`. Each group must be shorter than half the box. When unconnected: The "id" attribute if available, otherwise the index.
    place_groups : InputBoolean
        Wrap the centre of each group into the box. Otherwise each group stays nearest its average position
    wrap : InputMenu | Literal["Box", "Compact"]
        With Place Groups, wrap the centre of each group into the box spanned by the box vectors, or into the compact cell around the box centre, as Periodic Offset does
    update : InputBoolean
        Update the box lengths and angles with the simulation. As a modifier input, reads the "True" attribute by default.
    a : InputFloat
        a
    b : InputFloat
        b
    c : InputFloat
        c
    alpha : InputFloat
        alpha
    beta : InputFloat
        beta
    gamma : InputFloat
        gamma

    Inputs
    ------
    i.position : VectorSocket
        Positions to unwrap
    i.group_id : IntegerSocket
        Atoms with the same ID are kept whole, such as `fragment_id` for bonded molecules or `chain_id`. Each group must be shorter than half the box
    i.place_groups : BooleanSocket
        Wrap the centre of each group into the box. Otherwise each group stays nearest its average position
    i.wrap : MenuSocket
        With Place Groups, wrap the centre of each group into the box spanned by the box vectors, or into the compact cell around the box centre, as Periodic Offset does
    i.update : BooleanSocket
        Update the box lengths and angles with the simulation
    i.a : FloatSocket
        a
    i.b : FloatSocket
        b
    i.c : FloatSocket
        c
    i.alpha : FloatSocket
        alpha
    i.beta : FloatSocket
        beta
    i.gamma : FloatSocket
        gamma

    Outputs
    -------
    o.position : VectorSocket
        Positions with each group whole, unchanged without a box
    """

    _name = "Periodic Unwrap"
    _asset_name = "Periodic Unwrap"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Make groups of atoms whole across the periodic box, such as molecules split by its boundaries"
    }

    class _Inputs(SocketAccessor):
        position: VectorSocket
        """Positions to unwrap"""
        group_id: IntegerSocket
        """Atoms with the same ID are kept whole, such as `fragment_id` for bonded molecules or `chain_id`. Each group must be shorter than half the box"""
        place_groups: BooleanSocket
        """Wrap the centre of each group into the box. Otherwise each group stays nearest its average position"""
        wrap: MenuSocket
        """With Place Groups, wrap the centre of each group into the box spanned by the box vectors, or into the compact cell around the box centre, as Periodic Offset does"""
        update: BooleanSocket
        """Update the box lengths and angles with the simulation"""
        a: FloatSocket
        b: FloatSocket
        c: FloatSocket
        alpha: FloatSocket
        beta: FloatSocket
        gamma: FloatSocket

    class _Outputs(SocketAccessor):
        position: VectorSocket
        """Positions with each group whole, unchanged without a box"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        position: InputVector = Default.POSITION,
        group_id: InputInteger = Default.ID_OR_INDEX,
        place_groups: InputBoolean = False,
        wrap: InputMenu | Literal["Box", "Compact"] = "Box",
        update: InputBoolean = True,
        a: InputFloat = 0.0,
        b: InputFloat = 0.0,
        c: InputFloat = 0.0,
        alpha: InputFloat = 0.0,
        beta: InputFloat = 0.0,
        gamma: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Position": position,
                "Group ID": group_id,
                "Place Groups": place_groups,
                "Wrap": wrap,
                "Update": update,
                "a": a,
                "b": b,
                "c": c,
                "alpha": alpha,
                "beta": beta,
                "gamma": gamma,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        position = tree.inputs.vector(
            "Position",
            (0.0, 0.0, 0.0),
            description="Positions to unwrap",
            hide_value=True,
            default_input="POSITION",
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="Atoms with the same ID are kept whole, such as `fragment_id` for bonded molecules or `chain_id`. Each group must be shorter than half the box",
            hide_value=True,
            default_input="ID_OR_INDEX",
        )
        place_groups = tree.inputs.boolean(
            "Place Groups",
            False,
            description="Wrap the centre of each group into the box. Otherwise each group stays nearest its average position",
        )
        wrap = tree.inputs.menu(
            "Wrap",
            "Box",
            description="With Place Groups, wrap the centre of each group into the box spanned by the box vectors, or into the compact cell around the box centre, as Periodic Offset does",
            expanded=True,
        )
        with tree.inputs.panel("Periodic Box", default_closed=True):
            _update = tree.inputs.boolean(
                "Update",
                True,
                description="Update the box lengths and angles with the simulation",
                default_attribute="True",
            )
            with tree.inputs.panel("Lengths"):
                a = tree.inputs.float("a", 0.0)
                b = tree.inputs.float("b", 0.0)
                c_ = tree.inputs.float("c", 0.0)
            with tree.inputs.panel("Angles"):
                alpha = tree.inputs.float("alpha", 0.0)
                beta = tree.inputs.float("beta", 0.0)
                gamma = tree.inputs.float("gamma", 0.0)
        position_1 = tree.outputs.vector(
            "Position",
            description="Positions with each group whole, unchanged without a box",
        )

        periodic_box = PeriodicBox(a=a, b=b, c=c_, alpha=alpha, beta=beta, gamma=gamma)
        invert_matrix = g.InvertMatrix(periodic_box.o.matrix)
        transform_point = position.transform(invert_matrix.o.matrix)
        with g.Frame("Circular Mean"):
            _string = g.String(
                "The centre of each group in fractional coordinates, as the circular mean along each box vector, so a group split across a boundary is centred correctly."
            )
            vector_math = transform_point * math.tau
            mean = vector_math.cos().point.mean(group_id)
            mean_1 = vector_math.sin().point.mean(group_id)
            combine_xyz = g.CombineXYZ(
                x=mean_1.x.atan2(mean.x),
                y=mean_1.y.atan2(mean.y),
                z=mean_1.z.atan2(mean.z),
            )
            vector_math_1 = (combine_xyz.o.vector * 0.15915494).fraction()
        with g.Frame("Keep In Place"):
            _string_1 = g.String(
                "Without Place Groups, move the centre to its image nearest the group's plain average, so groups that are already whole stay where they are."
            )
            transform_point_1 = transform_point.point.mean(group_id).transform(
                periodic_box.o.matrix
            )
            transform_point_2 = vector_math_1.transform(periodic_box.o.matrix)
            periodic_nearest_image = PeriodicNearestImage(
                vector=transform_point_2 - transform_point_1,
                a=periodic_box.o.a,
                b=periodic_box.o.b,
                c=periodic_box.o.c,
            )
            vector_math_2 = transform_point_1 + periodic_nearest_image
        with g.Frame("Place Groups"):
            _string_2 = g.String(
                "With Place Groups, wrap the centre of each group into the box, or into the compact cell around the box centre."
            )
            periodic_offset = PeriodicOffset(
                position=transform_point_2,
                wrap=wrap,
                a=a,
                b=b,
                c=c_,
                alpha=alpha,
                beta=beta,
                gamma=gamma,
            )
        switch = place_groups.switch.vector(vector_math_2, periodic_offset)
        with g.Frame("Minimum Image"):
            _string_3 = g.String(
                "Compact wrapping the offset from the centre, moved to the middle of the box by half a box, then removing that half box, gives each atom's image nearest the centre, exactly for triclinic boxes too."
            )
            periodic_offset_1 = PeriodicOffset(
                position=position - switch,
                offset_factor=(0.5, 0.5, 0.5),
                wrap="Compact",
                a=a,
                b=b,
                c=c_,
                alpha=alpha,
                beta=beta,
                gamma=gamma,
            )
            vector_math_3 = periodic_offset_1.o.position - g.TransformDirection(
                transform=periodic_box.o.matrix, direction=(0.5, 0.5, 0.5)
            )
        (
            invert_matrix.o.invertible.switch.vector(position, switch + vector_math_3)
            >> position_1
        )


ASSET = PeriodicUnwrap

ASSET_METADATA = {
    "description": "Make groups of atoms whole across the periodic box, such as molecules split by its boundaries",
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
