# Node-group asset "Periodic Offset" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import Default, InputBoolean, InputFloat, InputMenu, InputVector
from ._shared.periodic_nearest_image import PeriodicNearestImage
from .periodic_box import PeriodicBox


class PeriodicOffset(AssetGeometryGroup):
    """
    Offset positions through the periodic box and wrap them back inside it

    Parameters
    ----------
    position : InputVector
        Positions to offset and wrap into the periodic box. When unconnected: The position from the context.
    offset_type : InputMenu | Literal["Factor", "Distance"]
        Offset by fractions of the box vectors, or by a distance
    offset_factor : InputVector
        Offset as a factor of the box vectors a, b and c, so 0.5 moves by half the box
    offset_distance : InputVector
        Offset as a distance in world units, where 1 m is 1 nm
    wrap : InputMenu | Literal["Box", "Compact"]
        Wrap into the box spanned by the box vectors, or into the compact cell around the box centre, where each atom is at its image nearest the centre (GROMACS `-ur compact`), such as a hexagonal prism for a hexagonal box
    update : InputBoolean
        Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node. As a modifier input, reads the "True" attribute by default. As a modifier input, reads the "True" attribute by default.
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
        Positions to offset and wrap into the periodic box
    i.offset_type : MenuSocket
        Offset by fractions of the box vectors, or by a distance
    i.offset_factor : VectorSocket
        Offset as a factor of the box vectors a, b and c, so 0.5 moves by half the box
    i.offset_distance : VectorSocket
        Offset as a distance in world units, where 1 m is 1 nm
    i.wrap : MenuSocket
        Wrap into the box spanned by the box vectors, or into the compact cell around the box centre, where each atom is at its image nearest the centre (GROMACS `-ur compact`), such as a hexagonal prism for a hexagonal box
    i.update : BooleanSocket
        Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node. As a modifier input, reads the "True" attribute by default.
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
        Offset positions wrapped into the periodic box, unchanged without a box
    """

    _name = "Periodic Offset"
    _asset_name = "Periodic Offset"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Offset positions through the periodic box and wrap them back inside it"
    }

    class _Inputs(SocketAccessor):
        position: VectorSocket
        """Positions to offset and wrap into the periodic box"""
        offset_type: MenuSocket
        """Offset by fractions of the box vectors, or by a distance"""
        offset_factor: VectorSocket
        """Offset as a factor of the box vectors a, b and c, so 0.5 moves by half the box"""
        offset_distance: VectorSocket
        """Offset as a distance in world units, where 1 m is 1 nm"""
        wrap: MenuSocket
        """Wrap into the box spanned by the box vectors, or into the compact cell around the box centre, where each atom is at its image nearest the centre (GROMACS `-ur compact`), such as a hexagonal prism for a hexagonal box"""
        update: BooleanSocket
        """Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node. As a modifier input, reads the "True" attribute by default."""
        a: FloatSocket
        b: FloatSocket
        c: FloatSocket
        alpha: FloatSocket
        beta: FloatSocket
        gamma: FloatSocket

    class _Outputs(SocketAccessor):
        position: VectorSocket
        """Offset positions wrapped into the periodic box, unchanged without a box"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        position: InputVector = Default.POSITION,
        offset_type: InputMenu | Literal["Factor", "Distance"] = "Factor",
        offset_factor: InputVector = None,
        offset_distance: InputVector = None,
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
                "Offset Type": offset_type,
                "Offset Factor": offset_factor,
                "Offset Distance": offset_distance,
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
            description="Positions to offset and wrap into the periodic box",
            hide_value=True,
            default_input="POSITION",
        )
        offset_type = tree.inputs.menu(
            "Offset Type",
            "Factor",
            description="Offset by fractions of the box vectors, or by a distance",
            expanded=True,
            optional_label=True,
        )
        offset_factor = tree.inputs.vector(
            "Offset Factor",
            (0.0, 0.0, 0.0),
            description="Offset as a factor of the box vectors a, b and c, so 0.5 moves by half the box",
            min_value=-1.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        offset_distance = tree.inputs.vector(
            "Offset Distance",
            (0.0, 0.0, 0.0),
            description="Offset as a distance in world units, where 1 m is 1 nm",
            subtype="TRANSLATION",
        )
        wrap = tree.inputs.menu(
            "Wrap",
            "Box",
            description="Wrap into the box spanned by the box vectors, or into the compact cell around the box centre, where each atom is at its image nearest the centre (GROMACS `-ur compact`), such as a hexagonal prism for a hexagonal box",
            expanded=True,
            optional_label=True,
        )
        with tree.inputs.panel("Periodic Box", default_closed=True):
            _update = tree.inputs.boolean(
                "Update",
                True,
                description='Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node. As a modifier input, reads the "True" attribute by default.',
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
            description="Offset positions wrapped into the periodic box, unchanged without a box",
        )

        periodic_box = PeriodicBox(a=a, b=b, c=c_, alpha=alpha, beta=beta, gamma=gamma)
        invert_matrix = g.InvertMatrix(periodic_box.o.matrix)
        menu_switch = g.MenuSwitch.vector(
            offset_type,
            {
                "Factor": position.transform(invert_matrix.o.matrix) + offset_factor,
                "Distance": (position + offset_distance).transform(
                    invert_matrix.o.matrix
                ),
            },
        )
        transform_point = menu_switch.o.output.fraction().transform(
            periodic_box.o.matrix
        )
        with g.Frame("Compact"):
            _string = g.String(
                "Move each position to its image nearest the centre of the box."
            )
            transform_point_1 = g.TransformPoint(
                transform=periodic_box.o.matrix, vector=(0.5, 0.5, 0.5)
            )
            periodic_nearest_image = PeriodicNearestImage(
                vector=transform_point - transform_point_1,
                a=periodic_box.o.a,
                b=periodic_box.o.b,
                c=periodic_box.o.c,
            )
            vector_math = transform_point_1.o.vector + periodic_nearest_image
        (
            invert_matrix.o.invertible.switch.vector(
                position,
                g.MenuSwitch.vector(
                    wrap, {"Box": transform_point, "Compact": vector_math}
                ).o.output,
            )
            >> position_1
        )


ASSET = PeriodicOffset

ASSET_METADATA = {
    "description": "Offset positions through the periodic box and wrap them back inside it",
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
