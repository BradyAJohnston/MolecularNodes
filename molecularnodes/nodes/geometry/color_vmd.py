# Node-group asset "Color VMD" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputMenu


class ColorVMD(AssetGeometryGroup):
    """
    Color VMD

    Parameters
    ----------
    color : InputMenu | Literal["blue", "red", "gray", "orange", "yellow", "tan", "silver", "green", "white", "pink", "cyan", "purple", "lime", "mauve", "ochre", "iceblue", "black", "yellow2", "yellow3", "green2", "green3", "cyan2", "cyan3", "blue2", "blue3", "violet", "violet2", "magenta", "magenta2", "red2", "red3", "orange2", "orange3"]
        VMD colour, in colour ID order (0 to 32)

    Inputs
    ------
    i.color : MenuSocket
        VMD colour, in colour ID order (0 to 32)

    Outputs
    -------
    o.color : ColorSocket
        The named colour, in linear RGB
    """

    _name = "Color VMD"
    _asset_name = "Color VMD"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        color: MenuSocket
        """VMD colour, in colour ID order (0 to 32)"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The named colour, in linear RGB"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        color: InputMenu
        | Literal[
            "blue",
            "red",
            "gray",
            "orange",
            "yellow",
            "tan",
            "silver",
            "green",
            "white",
            "pink",
            "cyan",
            "purple",
            "lime",
            "mauve",
            "ochre",
            "iceblue",
            "black",
            "yellow2",
            "yellow3",
            "green2",
            "green3",
            "cyan2",
            "cyan3",
            "blue2",
            "blue3",
            "violet",
            "violet2",
            "magenta",
            "magenta2",
            "red2",
            "red3",
            "orange2",
            "orange3",
        ] = "blue",
    ):
        super().__init__(**{"Color": color})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        color = tree.inputs.menu(
            "Color", "blue", description="VMD colour, in colour ID order (0 to 32)"
        )
        color_1 = tree.outputs.color(
            "Color", (1.0, 1.0, 1.0, 1.0), description="The named colour, in linear RGB"
        )

        (
            g.MenuSwitch.color(
                color,
                {
                    "blue": ((0.0, 0.0, 1.0, 1.0), "VMD blue, #0000FF"),
                    "red": ((1.0, 0.0, 0.0, 1.0), "VMD red, #FF0000"),
                    "gray": ((0.100482, 0.100482, 0.100482, 1.0), "VMD gray, #595959"),
                    "orange": ((1.0, 0.214041, 0.0, 1.0), "VMD orange, #FF8000"),
                    "yellow": ((1.0, 1.0, 0.0, 1.0), "VMD yellow, #FFFF00"),
                    "tan": ((0.214041, 0.214041, 0.033105, 1.0), "VMD tan, #808033"),
                    "silver": (
                        (0.318547, 0.318547, 0.318547, 1.0),
                        "VMD silver, #999999",
                    ),
                    "green": ((0.0, 1.0, 0.0, 1.0), "VMD green, #00FF00"),
                    "white": ((1.0, 1.0, 1.0, 1.0), "VMD white, #FFFFFF"),
                    "pink": ((1.0, 0.318547, 0.318547, 1.0), "VMD pink, #FF9999"),
                    "cyan": ((0.050876, 0.522522, 0.522522, 1.0), "VMD cyan, #40BFBF"),
                    "purple": ((0.380056, 0.0, 0.380056, 1.0), "VMD purple, #A600A6"),
                    "lime": ((0.214041, 0.787412, 0.132868, 1.0), "VMD lime, #80E666"),
                    "mauve": (
                        (0.787412, 0.132868, 0.447988, 1.0),
                        "VMD mauve, #E666B2",
                    ),
                    "ochre": ((0.214041, 0.073239, 0.0, 1.0), "VMD ochre, #804C00"),
                    "iceblue": (
                        (0.214041, 0.214041, 0.522522, 1.0),
                        "VMD iceblue, #8080BF",
                    ),
                    "black": ((0.0, 0.0, 0.0, 1.0), "VMD black, #000000"),
                    "yellow2": (
                        (0.748414, 0.933107, 0.001548, 1.0),
                        "VMD yellow2, #E0F705",
                    ),
                    "yellow3": (
                        (0.263273, 0.787412, 0.001548, 1.0),
                        "VMD yellow3, #8CE605",
                    ),
                    "green2": ((0.0, 0.787412, 0.003096, 1.0), "VMD green2, #00E60A"),
                    "green3": ((0.0, 0.787412, 0.214041, 1.0), "VMD green3, #00E680"),
                    "cyan2": ((0.0, 0.748414, 1.0, 1.0), "VMD cyan2, #00E0FF"),
                    "cyan3": ((0.0, 0.538236, 1.0, 1.0), "VMD cyan3, #00C2FF"),
                    "blue2": ((0.001548, 0.11928, 0.406448, 1.0), "VMD blue2, #0561AB"),
                    "blue3": (
                        (0.000774, 0.003096, 0.848088, 1.0),
                        "VMD blue3, #030AED",
                    ),
                    "violet": ((0.059254, 0.0, 0.955105, 1.0), "VMD violet, #4500FA"),
                    "violet2": ((0.170645, 0.0, 0.787412, 1.0), "VMD violet2, #7300E6"),
                    "magenta": ((0.787412, 0.0, 0.787412, 1.0), "VMD magenta, #E600E6"),
                    "magenta2": ((1.0, 0.0, 0.393123, 1.0), "VMD magenta2, #FF00A8"),
                    "red2": ((0.955105, 0.0, 0.043234, 1.0), "VMD red2, #FA003B"),
                    "red3": ((0.620916, 0.0, 0.0, 1.0), "VMD red3, #CF0000"),
                    "orange2": ((0.767769, 0.100482, 0.0, 1.0), "VMD orange2, #E35900"),
                    "orange3": ((0.911408, 0.477, 0.0, 1.0), "VMD orange3, #F5B800"),
                },
            )
            >> color_1
        )


ASSET = ColorVMD

ASSET_METADATA = {
    "description": "A named colour from VMD",
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
