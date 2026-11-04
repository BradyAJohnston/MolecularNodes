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
        Positions to move and wrap into the periodic box. When unconnected: The position from the context.
    offset_type : InputMenu | Literal["Factor", "Distance"]
        Give the offset as a factor of the box vectors, or as a distance
    offset_factor : InputVector
        Offset as a factor of box vectors A, B and C, so 0.5 moves by half the box along each
    offset_distance : InputVector
        Offset as a distance in world units, where 1 m is 1 nm
    wrap : InputMenu | Literal["Box", "Compact"]
        Wrap into the box spanned by the box vectors, or into the compact cell around the box centre, putting each position at its image nearest the centre (GROMACS `-ur compact`). For a hexagonal box the compact cell is a hexagonal prism
    update : InputBoolean
        Set the box lengths and angles from the trajectory on every frame change, as the molecule updates its positions. Turn off to set them by hand. As a modifier input, reads the "True" attribute by default.
    a : InputFloat
        Length of box vector A in world units, where 1 m is 1 nm
    b : InputFloat
        Length of box vector B in world units, where 1 m is 1 nm
    c : InputFloat
        Length of box vector C in world units, where 1 m is 1 nm
    alpha : InputFloat
        Angle between box vectors B and C, in degrees
    beta : InputFloat
        Angle between box vectors A and C, in degrees
    gamma : InputFloat
        Angle between box vectors A and B, in degrees

    Inputs
    ------
    i.position : VectorSocket
        Positions to move and wrap into the periodic box
    i.offset_type : MenuSocket
        Give the offset as a factor of the box vectors, or as a distance
    i.offset_factor : VectorSocket
        Offset as a factor of box vectors A, B and C, so 0.5 moves by half the box along each
    i.offset_distance : VectorSocket
        Offset as a distance in world units, where 1 m is 1 nm
    i.wrap : MenuSocket
        Wrap into the box spanned by the box vectors, or into the compact cell around the box centre, putting each position at its image nearest the centre (GROMACS `-ur compact`). For a hexagonal box the compact cell is a hexagonal prism
    i.update : BooleanSocket
        Set the box lengths and angles from the trajectory on every frame change, as the molecule updates its positions. Turn off to set them by hand
    i.a : FloatSocket
        Length of box vector A in world units, where 1 m is 1 nm
    i.b : FloatSocket
        Length of box vector B in world units, where 1 m is 1 nm
    i.c : FloatSocket
        Length of box vector C in world units, where 1 m is 1 nm
    i.alpha : FloatSocket
        Angle between box vectors B and C, in degrees
    i.beta : FloatSocket
        Angle between box vectors A and C, in degrees
    i.gamma : FloatSocket
        Angle between box vectors A and B, in degrees

    Outputs
    -------
    o.position : VectorSocket
        Positions offset and wrapped into the periodic box, or unchanged without a box
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
        """Positions to move and wrap into the periodic box"""
        offset_type: MenuSocket
        """Give the offset as a factor of the box vectors, or as a distance"""
        offset_factor: VectorSocket
        """Offset as a factor of box vectors A, B and C, so 0.5 moves by half the box along each"""
        offset_distance: VectorSocket
        """Offset as a distance in world units, where 1 m is 1 nm"""
        wrap: MenuSocket
        """Wrap into the box spanned by the box vectors, or into the compact cell around the box centre, putting each position at its image nearest the centre (GROMACS `-ur compact`). For a hexagonal box the compact cell is a hexagonal prism"""
        update: BooleanSocket
        """Set the box lengths and angles from the trajectory on every frame change, as the molecule updates its positions. Turn off to set them by hand"""
        a: FloatSocket
        """Length of box vector A in world units, where 1 m is 1 nm"""
        b: FloatSocket
        """Length of box vector B in world units, where 1 m is 1 nm"""
        c: FloatSocket
        """Length of box vector C in world units, where 1 m is 1 nm"""
        alpha: FloatSocket
        """Angle between box vectors B and C, in degrees"""
        beta: FloatSocket
        """Angle between box vectors A and C, in degrees"""
        gamma: FloatSocket
        """Angle between box vectors A and B, in degrees"""

    class _Outputs(SocketAccessor):
        position: VectorSocket
        """Positions offset and wrapped into the periodic box, or unchanged without a box"""

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
            description="Positions to move and wrap into the periodic box",
            hide_value=True,
            default_input="POSITION",
        )
        offset_type = tree.inputs.menu(
            "Offset Type",
            "Factor",
            description="Give the offset as a factor of the box vectors, or as a distance",
            expanded=True,
            optional_label=True,
        )
        offset_factor = tree.inputs.vector(
            "Offset Factor",
            (0.0, 0.0, 0.0),
            description="Offset as a factor of box vectors A, B and C, so 0.5 moves by half the box along each",
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
            description="Wrap into the box spanned by the box vectors, or into the compact cell around the box centre, putting each position at its image nearest the centre (GROMACS `-ur compact`). For a hexagonal box the compact cell is a hexagonal prism",
            expanded=True,
            optional_label=True,
        )
        with tree.inputs.panel("Periodic Box", default_closed=True):
            _update = tree.inputs.boolean(
                "Update",
                True,
                description="Set the box lengths and angles from the trajectory on every frame change, as the molecule updates its positions. Turn off to set them by hand",
                default_attribute="True",
            )
            with tree.inputs.panel("Lengths"):
                a = tree.inputs.float(
                    "a",
                    0.0,
                    description="Length of box vector A in world units, where 1 m is 1 nm",
                    subtype="DISTANCE",
                )
                b = tree.inputs.float(
                    "b",
                    0.0,
                    description="Length of box vector B in world units, where 1 m is 1 nm",
                    subtype="DISTANCE",
                )
                c_ = tree.inputs.float(
                    "c",
                    0.0,
                    description="Length of box vector C in world units, where 1 m is 1 nm",
                    subtype="DISTANCE",
                )
            with tree.inputs.panel("Angles"):
                alpha = tree.inputs.float(
                    "alpha",
                    0.0,
                    description="Angle between box vectors B and C, in degrees",
                )
                beta = tree.inputs.float(
                    "beta",
                    0.0,
                    description="Angle between box vectors A and C, in degrees",
                )
                gamma = tree.inputs.float(
                    "gamma",
                    0.0,
                    description="Angle between box vectors A and B, in degrees",
                )
        position_1 = tree.outputs.vector(
            "Position",
            description="Positions offset and wrapped into the periodic box, or unchanged without a box",
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
