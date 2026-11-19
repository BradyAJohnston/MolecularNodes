# Node-group asset "Color CSS" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class ColorCSS(AssetGeometryGroup):
    """
    Color CSS

    Parameters
    ----------
    color : InputMenu | Literal["aliceblue", "antiquewhite", "aqua", "aquamarine", "azure", "beige", "bisque", "black", "blanchedalmond", "blue", "blueviolet", "brown", "burlywood", "cadetblue", "chartreuse", "chocolate", "coral", "cornflowerblue", "cornsilk", "crimson", "cyan", "darkblue", "darkcyan", "darkgoldenrod", "darkgray", "darkgreen", "darkgrey", "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange", "darkorchid", "darkred", "darksalmon", "darkseagreen", "darkslateblue", "darkslategray", "darkslategrey", "darkturquoise", "darkviolet", "deeppink", "deepskyblue", "dimgray", "dimgrey", "dodgerblue", "firebrick", "floralwhite", "forestgreen", "fuchsia", "gainsboro", "ghostwhite", "gold", "goldenrod", "gray", "green", "greenyellow", "grey", "honeydew", "hotpink", "indianred", "indigo", "ivory", "khaki", "lavender", "lavenderblush", "lawngreen", "lemonchiffon", "lightblue", "lightcoral", "lightcyan", "lightgoldenrodyellow", "lightgray", "lightgreen", "lightgrey", "lightpink", "lightsalmon", "lightseagreen", "lightskyblue", "lightslategray", "lightslategrey", "lightsteelblue", "lightyellow", "lime", "limegreen", "linen", "magenta", "maroon", "mediumaquamarine", "mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen", "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue", "mintcream", "mistyrose", "moccasin", "navajowhite", "navy", "oldlace", "olive", "olivedrab", "orange", "orangered", "orchid", "palegoldenrod", "palegreen", "paleturquoise", "palevioletred", "papayawhip", "peachpuff", "peru", "pink", "plum", "powderblue", "purple", "rebeccapurple", "red", "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown", "seagreen", "seashell", "sienna", "silver", "skyblue", "slateblue", "slategray", "slategrey", "snow", "springgreen", "steelblue", "tan", "teal", "thistle", "tomato", "turquoise", "violet", "wheat", "white", "whitesmoke", "yellow", "yellowgreen"]
        CSS colour name, as used by MolViewSpec and ChimeraX

    Inputs
    ------
    i.color : MenuSocket
        CSS colour name, as used by MolViewSpec and ChimeraX

    Outputs
    -------
    o.color : ColorSocket
        The named colour, in linear RGB
    """

    _name = "Color CSS"
    _asset_name = "Color CSS"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        color: MenuSocket
        """CSS colour name, as used by MolViewSpec and ChimeraX"""

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
            "aliceblue",
            "antiquewhite",
            "aqua",
            "aquamarine",
            "azure",
            "beige",
            "bisque",
            "black",
            "blanchedalmond",
            "blue",
            "blueviolet",
            "brown",
            "burlywood",
            "cadetblue",
            "chartreuse",
            "chocolate",
            "coral",
            "cornflowerblue",
            "cornsilk",
            "crimson",
            "cyan",
            "darkblue",
            "darkcyan",
            "darkgoldenrod",
            "darkgray",
            "darkgreen",
            "darkgrey",
            "darkkhaki",
            "darkmagenta",
            "darkolivegreen",
            "darkorange",
            "darkorchid",
            "darkred",
            "darksalmon",
            "darkseagreen",
            "darkslateblue",
            "darkslategray",
            "darkslategrey",
            "darkturquoise",
            "darkviolet",
            "deeppink",
            "deepskyblue",
            "dimgray",
            "dimgrey",
            "dodgerblue",
            "firebrick",
            "floralwhite",
            "forestgreen",
            "fuchsia",
            "gainsboro",
            "ghostwhite",
            "gold",
            "goldenrod",
            "gray",
            "green",
            "greenyellow",
            "grey",
            "honeydew",
            "hotpink",
            "indianred",
            "indigo",
            "ivory",
            "khaki",
            "lavender",
            "lavenderblush",
            "lawngreen",
            "lemonchiffon",
            "lightblue",
            "lightcoral",
            "lightcyan",
            "lightgoldenrodyellow",
            "lightgray",
            "lightgreen",
            "lightgrey",
            "lightpink",
            "lightsalmon",
            "lightseagreen",
            "lightskyblue",
            "lightslategray",
            "lightslategrey",
            "lightsteelblue",
            "lightyellow",
            "lime",
            "limegreen",
            "linen",
            "magenta",
            "maroon",
            "mediumaquamarine",
            "mediumblue",
            "mediumorchid",
            "mediumpurple",
            "mediumseagreen",
            "mediumslateblue",
            "mediumspringgreen",
            "mediumturquoise",
            "mediumvioletred",
            "midnightblue",
            "mintcream",
            "mistyrose",
            "moccasin",
            "navajowhite",
            "navy",
            "oldlace",
            "olive",
            "olivedrab",
            "orange",
            "orangered",
            "orchid",
            "palegoldenrod",
            "palegreen",
            "paleturquoise",
            "palevioletred",
            "papayawhip",
            "peachpuff",
            "peru",
            "pink",
            "plum",
            "powderblue",
            "purple",
            "rebeccapurple",
            "red",
            "rosybrown",
            "royalblue",
            "saddlebrown",
            "salmon",
            "sandybrown",
            "seagreen",
            "seashell",
            "sienna",
            "silver",
            "skyblue",
            "slateblue",
            "slategray",
            "slategrey",
            "snow",
            "springgreen",
            "steelblue",
            "tan",
            "teal",
            "thistle",
            "tomato",
            "turquoise",
            "violet",
            "wheat",
            "white",
            "whitesmoke",
            "yellow",
            "yellowgreen",
        ] = "steelblue",
    ):
        super().__init__(**{"Color": color})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        color = tree.inputs.menu(
            "Color",
            "steelblue",
            description="CSS colour name, as used by MolViewSpec and ChimeraX",
        )
        color_1 = tree.outputs.color(
            "Color", (1.0, 1.0, 1.0, 1.0), description="The named colour, in linear RGB"
        )

        (
            g.MenuSwitch.color(
                color,
                {
                    "aliceblue": (
                        (0.871367, 0.938686, 1.0, 1.0),
                        "CSS aliceblue, #F0F8FF",
                    ),
                    "antiquewhite": (
                        (0.955973, 0.83077, 0.679542, 1.0),
                        "CSS antiquewhite, #FAEBD7",
                    ),
                    "aqua": ((0.0, 1.0, 1.0, 1.0), "CSS aqua, #00FFFF"),
                    "aquamarine": (
                        (0.212231, 1.0, 0.658375, 1.0),
                        "CSS aquamarine, #7FFFD4",
                    ),
                    "azure": ((0.871367, 1.0, 1.0, 1.0), "CSS azure, #F0FFFF"),
                    "beige": (
                        (0.913099, 0.913099, 0.715694, 1.0),
                        "CSS beige, #F5F5DC",
                    ),
                    "bisque": ((1.0, 0.775822, 0.552011, 1.0), "CSS bisque, #FFE4C4"),
                    "black": ((0.0, 0.0, 0.0, 1.0), "CSS black, #000000"),
                    "blanchedalmond": (
                        (1.0, 0.83077, 0.610496, 1.0),
                        "CSS blanchedalmond, #FFEBCD",
                    ),
                    "blue": ((0.0, 0.0, 1.0, 1.0), "CSS blue, #0000FF"),
                    "blueviolet": (
                        (0.254152, 0.024158, 0.760525, 1.0),
                        "CSS blueviolet, #8A2BE2",
                    ),
                    "brown": (
                        (0.376262, 0.023153, 0.023153, 1.0),
                        "CSS brown, #A52A2A",
                    ),
                    "burlywood": (
                        (0.730461, 0.47932, 0.242281, 1.0),
                        "CSS burlywood, #DEB887",
                    ),
                    "cadetblue": (
                        (0.114435, 0.341914, 0.351533, 1.0),
                        "CSS cadetblue, #5F9EA0",
                    ),
                    "chartreuse": (
                        (0.212231, 1.0, 0.0, 1.0),
                        "CSS chartreuse, #7FFF00",
                    ),
                    "chocolate": (
                        (0.64448, 0.141263, 0.012983, 1.0),
                        "CSS chocolate, #D2691E",
                    ),
                    "coral": ((1.0, 0.212231, 0.08022, 1.0), "CSS coral, #FF7F50"),
                    "cornflowerblue": (
                        (0.127438, 0.300544, 0.846873, 1.0),
                        "CSS cornflowerblue, #6495ED",
                    ),
                    "cornsilk": (
                        (1.0, 0.938686, 0.715694, 1.0),
                        "CSS cornsilk, #FFF8DC",
                    ),
                    "crimson": (
                        (0.715694, 0.006995, 0.045186, 1.0),
                        "CSS crimson, #DC143C",
                    ),
                    "cyan": ((0.0, 1.0, 1.0, 1.0), "CSS cyan, #00FFFF"),
                    "darkblue": ((0.0, 0.0, 0.258183, 1.0), "CSS darkblue, #00008B"),
                    "darkcyan": (
                        (0.0, 0.258183, 0.258183, 1.0),
                        "CSS darkcyan, #008B8B",
                    ),
                    "darkgoldenrod": (
                        (0.47932, 0.238398, 0.003347, 1.0),
                        "CSS darkgoldenrod, #B8860B",
                    ),
                    "darkgray": (
                        (0.396755, 0.396755, 0.396755, 1.0),
                        "CSS darkgray, #A9A9A9",
                    ),
                    "darkgreen": ((0.0, 0.127438, 0.0, 1.0), "CSS darkgreen, #006400"),
                    "darkgrey": (
                        (0.396755, 0.396755, 0.396755, 1.0),
                        "CSS darkgrey, #A9A9A9",
                    ),
                    "darkkhaki": (
                        (0.508881, 0.473531, 0.147027, 1.0),
                        "CSS darkkhaki, #BDB76B",
                    ),
                    "darkmagenta": (
                        (0.258183, 0.0, 0.258183, 1.0),
                        "CSS darkmagenta, #8B008B",
                    ),
                    "darkolivegreen": (
                        (0.090842, 0.147027, 0.028426, 1.0),
                        "CSS darkolivegreen, #556B2F",
                    ),
                    "darkorange": (
                        (1.0, 0.262251, 0.0, 1.0),
                        "CSS darkorange, #FF8C00",
                    ),
                    "darkorchid": (
                        (0.318547, 0.031896, 0.603827, 1.0),
                        "CSS darkorchid, #9932CC",
                    ),
                    "darkred": ((0.258183, 0.0, 0.0, 1.0), "CSS darkred, #8B0000"),
                    "darksalmon": (
                        (0.814847, 0.304987, 0.194618, 1.0),
                        "CSS darksalmon, #E9967A",
                    ),
                    "darkseagreen": (
                        (0.274677, 0.502886, 0.274677, 1.0),
                        "CSS darkseagreen, #8FBC8F",
                    ),
                    "darkslateblue": (
                        (0.064803, 0.046665, 0.258183, 1.0),
                        "CSS darkslateblue, #483D8B",
                    ),
                    "darkslategray": (
                        (0.028426, 0.078187, 0.078187, 1.0),
                        "CSS darkslategray, #2F4F4F",
                    ),
                    "darkslategrey": (
                        (0.028426, 0.078187, 0.078187, 1.0),
                        "CSS darkslategrey, #2F4F4F",
                    ),
                    "darkturquoise": (
                        (0.0, 0.617207, 0.637597, 1.0),
                        "CSS darkturquoise, #00CED1",
                    ),
                    "darkviolet": (
                        (0.296138, 0.0, 0.651406, 1.0),
                        "CSS darkviolet, #9400D3",
                    ),
                    "deeppink": (
                        (1.0, 0.006995, 0.291771, 1.0),
                        "CSS deeppink, #FF1493",
                    ),
                    "deepskyblue": (
                        (0.0, 0.520996, 1.0, 1.0),
                        "CSS deepskyblue, #00BFFF",
                    ),
                    "dimgray": (
                        (0.141263, 0.141263, 0.141263, 1.0),
                        "CSS dimgray, #696969",
                    ),
                    "dimgrey": (
                        (0.141263, 0.141263, 0.141263, 1.0),
                        "CSS dimgrey, #696969",
                    ),
                    "dodgerblue": (
                        (0.012983, 0.278894, 1.0, 1.0),
                        "CSS dodgerblue, #1E90FF",
                    ),
                    "firebrick": (
                        (0.445201, 0.015996, 0.015996, 1.0),
                        "CSS firebrick, #B22222",
                    ),
                    "floralwhite": (
                        (1.0, 0.955973, 0.871367, 1.0),
                        "CSS floralwhite, #FFFAF0",
                    ),
                    "forestgreen": (
                        (0.015996, 0.258183, 0.015996, 1.0),
                        "CSS forestgreen, #228B22",
                    ),
                    "fuchsia": ((1.0, 0.0, 1.0, 1.0), "CSS fuchsia, #FF00FF"),
                    "gainsboro": (
                        (0.715694, 0.715694, 0.715694, 1.0),
                        "CSS gainsboro, #DCDCDC",
                    ),
                    "ghostwhite": (
                        (0.938686, 0.938686, 1.0, 1.0),
                        "CSS ghostwhite, #F8F8FF",
                    ),
                    "gold": ((1.0, 0.679542, 0.0, 1.0), "CSS gold, #FFD700"),
                    "goldenrod": (
                        (0.701102, 0.376262, 0.014444, 1.0),
                        "CSS goldenrod, #DAA520",
                    ),
                    "gray": ((0.215861, 0.215861, 0.215861, 1.0), "CSS gray, #808080"),
                    "green": ((0.0, 0.215861, 0.0, 1.0), "CSS green, #008000"),
                    "greenyellow": (
                        (0.417885, 1.0, 0.028426, 1.0),
                        "CSS greenyellow, #ADFF2F",
                    ),
                    "grey": ((0.215861, 0.215861, 0.215861, 1.0), "CSS grey, #808080"),
                    "honeydew": (
                        (0.871367, 1.0, 0.871367, 1.0),
                        "CSS honeydew, #F0FFF0",
                    ),
                    "hotpink": ((1.0, 0.141263, 0.456411, 1.0), "CSS hotpink, #FF69B4"),
                    "indianred": (
                        (0.610496, 0.107023, 0.107023, 1.0),
                        "CSS indianred, #CD5C5C",
                    ),
                    "indigo": ((0.07036, 0.0, 0.223228, 1.0), "CSS indigo, #4B0082"),
                    "ivory": ((1.0, 1.0, 0.871367, 1.0), "CSS ivory, #FFFFF0"),
                    "khaki": (
                        (0.871367, 0.791298, 0.262251, 1.0),
                        "CSS khaki, #F0E68C",
                    ),
                    "lavender": (
                        (0.791298, 0.791298, 0.955973, 1.0),
                        "CSS lavender, #E6E6FA",
                    ),
                    "lavenderblush": (
                        (1.0, 0.871367, 0.913099, 1.0),
                        "CSS lavenderblush, #FFF0F5",
                    ),
                    "lawngreen": (
                        (0.201556, 0.973445, 0.0, 1.0),
                        "CSS lawngreen, #7CFC00",
                    ),
                    "lemonchiffon": (
                        (1.0, 0.955973, 0.610496, 1.0),
                        "CSS lemonchiffon, #FFFACD",
                    ),
                    "lightblue": (
                        (0.417885, 0.686685, 0.791298, 1.0),
                        "CSS lightblue, #ADD8E6",
                    ),
                    "lightcoral": (
                        (0.871367, 0.215861, 0.215861, 1.0),
                        "CSS lightcoral, #F08080",
                    ),
                    "lightcyan": ((0.745404, 1.0, 1.0, 1.0), "CSS lightcyan, #E0FFFF"),
                    "lightgoldenrodyellow": (
                        (0.955973, 0.955973, 0.64448, 1.0),
                        "CSS lightgoldenrodyellow, #FAFAD2",
                    ),
                    "lightgray": (
                        (0.651406, 0.651406, 0.651406, 1.0),
                        "CSS lightgray, #D3D3D3",
                    ),
                    "lightgreen": (
                        (0.278894, 0.854993, 0.278894, 1.0),
                        "CSS lightgreen, #90EE90",
                    ),
                    "lightgrey": (
                        (0.651406, 0.651406, 0.651406, 1.0),
                        "CSS lightgrey, #D3D3D3",
                    ),
                    "lightpink": (
                        (1.0, 0.467784, 0.533276, 1.0),
                        "CSS lightpink, #FFB6C1",
                    ),
                    "lightsalmon": (
                        (1.0, 0.351533, 0.194618, 1.0),
                        "CSS lightsalmon, #FFA07A",
                    ),
                    "lightseagreen": (
                        (0.014444, 0.445201, 0.401978, 1.0),
                        "CSS lightseagreen, #20B2AA",
                    ),
                    "lightskyblue": (
                        (0.242281, 0.617207, 0.955973, 1.0),
                        "CSS lightskyblue, #87CEFA",
                    ),
                    "lightslategray": (
                        (0.184475, 0.246201, 0.318547, 1.0),
                        "CSS lightslategray, #778899",
                    ),
                    "lightslategrey": (
                        (0.184475, 0.246201, 0.318547, 1.0),
                        "CSS lightslategrey, #778899",
                    ),
                    "lightsteelblue": (
                        (0.434154, 0.552011, 0.730461, 1.0),
                        "CSS lightsteelblue, #B0C4DE",
                    ),
                    "lightyellow": (
                        (1.0, 1.0, 0.745404, 1.0),
                        "CSS lightyellow, #FFFFE0",
                    ),
                    "lime": ((0.0, 1.0, 0.0, 1.0), "CSS lime, #00FF00"),
                    "limegreen": (
                        (0.031896, 0.610496, 0.031896, 1.0),
                        "CSS limegreen, #32CD32",
                    ),
                    "linen": (
                        (0.955973, 0.871367, 0.791298, 1.0),
                        "CSS linen, #FAF0E6",
                    ),
                    "magenta": ((1.0, 0.0, 1.0, 1.0), "CSS magenta, #FF00FF"),
                    "maroon": ((0.215861, 0.0, 0.0, 1.0), "CSS maroon, #800000"),
                    "mediumaquamarine": (
                        (0.132868, 0.610496, 0.401978, 1.0),
                        "CSS mediumaquamarine, #66CDAA",
                    ),
                    "mediumblue": (
                        (0.0, 0.0, 0.610496, 1.0),
                        "CSS mediumblue, #0000CD",
                    ),
                    "mediumorchid": (
                        (0.491021, 0.090842, 0.651406, 1.0),
                        "CSS mediumorchid, #BA55D3",
                    ),
                    "mediumpurple": (
                        (0.291771, 0.162029, 0.708376, 1.0),
                        "CSS mediumpurple, #9370DB",
                    ),
                    "mediumseagreen": (
                        (0.045186, 0.450786, 0.165132, 1.0),
                        "CSS mediumseagreen, #3CB371",
                    ),
                    "mediumslateblue": (
                        (0.198069, 0.138432, 0.854993, 1.0),
                        "CSS mediumslateblue, #7B68EE",
                    ),
                    "mediumspringgreen": (
                        (0.0, 0.955973, 0.323143, 1.0),
                        "CSS mediumspringgreen, #00FA9A",
                    ),
                    "mediumturquoise": (
                        (0.064803, 0.637597, 0.603827, 1.0),
                        "CSS mediumturquoise, #48D1CC",
                    ),
                    "mediumvioletred": (
                        (0.571125, 0.007499, 0.234551, 1.0),
                        "CSS mediumvioletred, #C71585",
                    ),
                    "midnightblue": (
                        (0.009721, 0.009721, 0.162029, 1.0),
                        "CSS midnightblue, #191970",
                    ),
                    "mintcream": (
                        (0.913099, 1.0, 0.955973, 1.0),
                        "CSS mintcream, #F5FFFA",
                    ),
                    "mistyrose": (
                        (1.0, 0.775822, 0.752942, 1.0),
                        "CSS mistyrose, #FFE4E1",
                    ),
                    "moccasin": (
                        (1.0, 0.775822, 0.462077, 1.0),
                        "CSS moccasin, #FFE4B5",
                    ),
                    "navajowhite": (
                        (1.0, 0.730461, 0.417885, 1.0),
                        "CSS navajowhite, #FFDEAD",
                    ),
                    "navy": ((0.0, 0.0, 0.215861, 1.0), "CSS navy, #000080"),
                    "oldlace": (
                        (0.982251, 0.913099, 0.791298, 1.0),
                        "CSS oldlace, #FDF5E6",
                    ),
                    "olive": ((0.215861, 0.215861, 0.0, 1.0), "CSS olive, #808000"),
                    "olivedrab": (
                        (0.147027, 0.270498, 0.016807, 1.0),
                        "CSS olivedrab, #6B8E23",
                    ),
                    "orange": ((1.0, 0.376262, 0.0, 1.0), "CSS orange, #FFA500"),
                    "orangered": ((1.0, 0.059511, 0.0, 1.0), "CSS orangered, #FF4500"),
                    "orchid": (
                        (0.701102, 0.162029, 0.672443, 1.0),
                        "CSS orchid, #DA70D6",
                    ),
                    "palegoldenrod": (
                        (0.854993, 0.806952, 0.401978, 1.0),
                        "CSS palegoldenrod, #EEE8AA",
                    ),
                    "palegreen": (
                        (0.313989, 0.964686, 0.313989, 1.0),
                        "CSS palegreen, #98FB98",
                    ),
                    "paleturquoise": (
                        (0.42869, 0.854993, 0.854993, 1.0),
                        "CSS paleturquoise, #AFEEEE",
                    ),
                    "palevioletred": (
                        (0.708376, 0.162029, 0.291771, 1.0),
                        "CSS palevioletred, #DB7093",
                    ),
                    "papayawhip": (
                        (1.0, 0.863157, 0.665387, 1.0),
                        "CSS papayawhip, #FFEFD5",
                    ),
                    "peachpuff": (
                        (1.0, 0.701102, 0.48515, 1.0),
                        "CSS peachpuff, #FFDAB9",
                    ),
                    "peru": ((0.610496, 0.234551, 0.049707, 1.0), "CSS peru, #CD853F"),
                    "pink": ((1.0, 0.527115, 0.597202, 1.0), "CSS pink, #FFC0CB"),
                    "plum": ((0.723055, 0.351533, 0.723055, 1.0), "CSS plum, #DDA0DD"),
                    "powderblue": (
                        (0.434154, 0.745404, 0.791298, 1.0),
                        "CSS powderblue, #B0E0E6",
                    ),
                    "purple": ((0.215861, 0.0, 0.215861, 1.0), "CSS purple, #800080"),
                    "rebeccapurple": (
                        (0.132868, 0.033105, 0.318547, 1.0),
                        "CSS rebeccapurple, #663399",
                    ),
                    "red": ((1.0, 0.0, 0.0, 1.0), "CSS red, #FF0000"),
                    "rosybrown": (
                        (0.502886, 0.274677, 0.274677, 1.0),
                        "CSS rosybrown, #BC8F8F",
                    ),
                    "royalblue": (
                        (0.052861, 0.141263, 0.752942, 1.0),
                        "CSS royalblue, #4169E1",
                    ),
                    "saddlebrown": (
                        (0.258183, 0.059511, 0.006512, 1.0),
                        "CSS saddlebrown, #8B4513",
                    ),
                    "salmon": (
                        (0.955973, 0.215861, 0.168269, 1.0),
                        "CSS salmon, #FA8072",
                    ),
                    "sandybrown": (
                        (0.904661, 0.371238, 0.116971, 1.0),
                        "CSS sandybrown, #F4A460",
                    ),
                    "seagreen": (
                        (0.027321, 0.258183, 0.095307, 1.0),
                        "CSS seagreen, #2E8B57",
                    ),
                    "seashell": (
                        (1.0, 0.913099, 0.854993, 1.0),
                        "CSS seashell, #FFF5EE",
                    ),
                    "sienna": (
                        (0.351533, 0.084376, 0.026241, 1.0),
                        "CSS sienna, #A0522D",
                    ),
                    "silver": (
                        (0.527115, 0.527115, 0.527115, 1.0),
                        "CSS silver, #C0C0C0",
                    ),
                    "skyblue": (
                        (0.242281, 0.617207, 0.83077, 1.0),
                        "CSS skyblue, #87CEEB",
                    ),
                    "slateblue": (
                        (0.144128, 0.102242, 0.610496, 1.0),
                        "CSS slateblue, #6A5ACD",
                    ),
                    "slategray": (
                        (0.162029, 0.215861, 0.278894, 1.0),
                        "CSS slategray, #708090",
                    ),
                    "slategrey": (
                        (0.162029, 0.215861, 0.278894, 1.0),
                        "CSS slategrey, #708090",
                    ),
                    "snow": ((1.0, 0.955973, 0.955973, 1.0), "CSS snow, #FFFAFA"),
                    "springgreen": (
                        (0.0, 1.0, 0.212231, 1.0),
                        "CSS springgreen, #00FF7F",
                    ),
                    "steelblue": (
                        (0.061246, 0.223228, 0.456411, 1.0),
                        "CSS steelblue, #4682B4",
                    ),
                    "tan": ((0.64448, 0.456411, 0.262251, 1.0), "CSS tan, #D2B48C"),
                    "teal": ((0.0, 0.215861, 0.215861, 1.0), "CSS teal, #008080"),
                    "thistle": (
                        (0.686685, 0.520996, 0.686685, 1.0),
                        "CSS thistle, #D8BFD8",
                    ),
                    "tomato": ((1.0, 0.124772, 0.06301, 1.0), "CSS tomato, #FF6347"),
                    "turquoise": (
                        (0.051269, 0.745404, 0.630757, 1.0),
                        "CSS turquoise, #40E0D0",
                    ),
                    "violet": (
                        (0.854993, 0.223228, 0.854993, 1.0),
                        "CSS violet, #EE82EE",
                    ),
                    "wheat": (
                        (0.913099, 0.730461, 0.450786, 1.0),
                        "CSS wheat, #F5DEB3",
                    ),
                    "white": ((1.0, 1.0, 1.0, 1.0), "CSS white, #FFFFFF"),
                    "whitesmoke": (
                        (0.913099, 0.913099, 0.913099, 1.0),
                        "CSS whitesmoke, #F5F5F5",
                    ),
                    "yellow": ((1.0, 1.0, 0.0, 1.0), "CSS yellow, #FFFF00"),
                    "yellowgreen": (
                        (0.323143, 0.610496, 0.031896, 1.0),
                        "CSS yellowgreen, #9ACD32",
                    ),
                },
            )
            >> color_1
        )


ASSET = ColorCSS

ASSET_METADATA = {
    "description": "A named CSS colour, as used by MolViewSpec and ChimeraX",
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
