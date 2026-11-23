# Node-group asset "Color PyMOL" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class ColorPyMOL(AssetGeometryGroup):
    """
    Color PyMOL

    Parameters
    ----------
    category : InputMenu | Literal["Reds", "Greens", "Blues", "Yellows", "Magentas", "Cyans", "Oranges", "Tints", "Grays", "Elements"]
        Group of colours, as in PyMOL's Color menu
    reds : InputMenu | Literal["red", "tv_red", "raspberry", "darksalmon", "salmon", "deepsalmon", "warmpink", "firebrick", "ruby", "chocolate", "brown"]
        PyMOL colour, used when Category is Reds
    greens : InputMenu | Literal["green", "tv_green", "chartreuse", "splitpea", "smudge", "palegreen", "limegreen", "lime", "limon", "forest"]
        PyMOL colour, used when Category is Greens
    blues : InputMenu | Literal["blue", "tv_blue", "marine", "slate", "lightblue", "skyblue", "purpleblue", "deepblue", "density"]
        PyMOL colour, used when Category is Blues
    yellows : InputMenu | Literal["yellow", "tv_yellow", "paleyellow", "yelloworange", "limon", "wheat", "sand"]
        PyMOL colour, used when Category is Yellows
    magentas : InputMenu | Literal["magenta", "lightmagenta", "hotpink", "pink", "lightpink", "dirtyviolet", "violet", "violetpurple", "purple", "deeppurple"]
        PyMOL colour, used when Category is Magentas
    cyans : InputMenu | Literal["cyan", "palecyan", "aquamarine", "greencyan", "teal", "deepteal", "lightteal"]
        PyMOL colour, used when Category is Cyans
    oranges : InputMenu | Literal["orange", "tv_orange", "brightorange", "lightorange", "yelloworange", "olive", "deepolive"]
        PyMOL colour, used when Category is Oranges
    tints : InputMenu | Literal["wheat", "palegreen", "lightblue", "paleyellow", "lightpink", "palecyan", "lightorange", "bluewhite"]
        PyMOL colour, used when Category is Tints
    grays : InputMenu | Literal["white", "gray90", "gray80", "gray70", "gray60", "gray50", "gray40", "gray30", "gray20", "gray10", "black"]
        PyMOL colour, used when Category is Grays
    elements : InputMenu | Literal["actinium", "aluminum", "americium", "antimony", "argon", "arsenic", "astatine", "barium", "berkelium", "beryllium", "bismuth", "bohrium", "boron", "bromine", "cadmium", "calcium", "californium", "carbon", "cerium", "cesium", "chlorine", "chromium", "cobalt", "copper", "curium", "deuterium", "dubnium", "dysprosium", "einsteinium", "erbium", "europium", "fermium", "fluorine", "francium", "gadolinium", "gallium", "germanium", "gold", "hafnium", "hassium", "helium", "holmium", "hydrogen", "indium", "iodine", "iridium", "iron", "krypton", "lanthanum", "lawrencium", "lead", "lithium", "lutetium", "magnesium", "manganese", "meitnerium", "mendelevium", "mercury", "molybdenum", "neodymium", "neon", "neptunium", "nickel", "niobium", "nitrogen", "nobelium", "osmium", "oxygen", "palladium", "phosphorus", "platinum", "plutonium", "polonium", "potassium", "praseodymium", "promethium", "protactinium", "radium", "radon", "rhenium", "rhodium", "rubidium", "ruthenium", "rutherfordium", "samarium", "scandium", "seaborgium", "selenium", "silicon", "silver", "sodium", "strontium", "sulfur", "tantalum", "technetium", "tellurium", "terbium", "thallium", "thorium", "thulium", "tin", "titanium", "tungsten", "uranium", "vanadium", "xenon", "ytterbium", "yttrium", "zinc", "zirconium"]
        PyMOL colour, used when Category is Elements

    Inputs
    ------
    i.category : MenuSocket
        Group of colours, as in PyMOL's Color menu
    i.reds : MenuSocket
        PyMOL colour, used when Category is Reds
    i.greens : MenuSocket
        PyMOL colour, used when Category is Greens
    i.blues : MenuSocket
        PyMOL colour, used when Category is Blues
    i.yellows : MenuSocket
        PyMOL colour, used when Category is Yellows
    i.magentas : MenuSocket
        PyMOL colour, used when Category is Magentas
    i.cyans : MenuSocket
        PyMOL colour, used when Category is Cyans
    i.oranges : MenuSocket
        PyMOL colour, used when Category is Oranges
    i.tints : MenuSocket
        PyMOL colour, used when Category is Tints
    i.grays : MenuSocket
        PyMOL colour, used when Category is Grays
    i.elements : MenuSocket
        PyMOL colour, used when Category is Elements

    Outputs
    -------
    o.color : ColorSocket
        The named colour, in linear RGB
    """

    _name = "Color PyMOL"
    _asset_name = "Color PyMOL"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        category: MenuSocket
        """Group of colours, as in PyMOL's Color menu"""
        reds: MenuSocket
        """PyMOL colour, used when Category is Reds"""
        greens: MenuSocket
        """PyMOL colour, used when Category is Greens"""
        blues: MenuSocket
        """PyMOL colour, used when Category is Blues"""
        yellows: MenuSocket
        """PyMOL colour, used when Category is Yellows"""
        magentas: MenuSocket
        """PyMOL colour, used when Category is Magentas"""
        cyans: MenuSocket
        """PyMOL colour, used when Category is Cyans"""
        oranges: MenuSocket
        """PyMOL colour, used when Category is Oranges"""
        tints: MenuSocket
        """PyMOL colour, used when Category is Tints"""
        grays: MenuSocket
        """PyMOL colour, used when Category is Grays"""
        elements: MenuSocket
        """PyMOL colour, used when Category is Elements"""

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
        category: InputMenu
        | Literal[
            "Reds",
            "Greens",
            "Blues",
            "Yellows",
            "Magentas",
            "Cyans",
            "Oranges",
            "Tints",
            "Grays",
            "Elements",
        ] = "Reds",
        reds: InputMenu
        | Literal[
            "red",
            "tv_red",
            "raspberry",
            "darksalmon",
            "salmon",
            "deepsalmon",
            "warmpink",
            "firebrick",
            "ruby",
            "chocolate",
            "brown",
        ] = "red",
        greens: InputMenu
        | Literal[
            "green",
            "tv_green",
            "chartreuse",
            "splitpea",
            "smudge",
            "palegreen",
            "limegreen",
            "lime",
            "limon",
            "forest",
        ] = "green",
        blues: InputMenu
        | Literal[
            "blue",
            "tv_blue",
            "marine",
            "slate",
            "lightblue",
            "skyblue",
            "purpleblue",
            "deepblue",
            "density",
        ] = "blue",
        yellows: InputMenu
        | Literal[
            "yellow",
            "tv_yellow",
            "paleyellow",
            "yelloworange",
            "limon",
            "wheat",
            "sand",
        ] = "yellow",
        magentas: InputMenu
        | Literal[
            "magenta",
            "lightmagenta",
            "hotpink",
            "pink",
            "lightpink",
            "dirtyviolet",
            "violet",
            "violetpurple",
            "purple",
            "deeppurple",
        ] = "magenta",
        cyans: InputMenu
        | Literal[
            "cyan",
            "palecyan",
            "aquamarine",
            "greencyan",
            "teal",
            "deepteal",
            "lightteal",
        ] = "cyan",
        oranges: InputMenu
        | Literal[
            "orange",
            "tv_orange",
            "brightorange",
            "lightorange",
            "yelloworange",
            "olive",
            "deepolive",
        ] = "orange",
        tints: InputMenu
        | Literal[
            "wheat",
            "palegreen",
            "lightblue",
            "paleyellow",
            "lightpink",
            "palecyan",
            "lightorange",
            "bluewhite",
        ] = "wheat",
        grays: InputMenu
        | Literal[
            "white",
            "gray90",
            "gray80",
            "gray70",
            "gray60",
            "gray50",
            "gray40",
            "gray30",
            "gray20",
            "gray10",
            "black",
        ] = "white",
        elements: InputMenu
        | Literal[
            "actinium",
            "aluminum",
            "americium",
            "antimony",
            "argon",
            "arsenic",
            "astatine",
            "barium",
            "berkelium",
            "beryllium",
            "bismuth",
            "bohrium",
            "boron",
            "bromine",
            "cadmium",
            "calcium",
            "californium",
            "carbon",
            "cerium",
            "cesium",
            "chlorine",
            "chromium",
            "cobalt",
            "copper",
            "curium",
            "deuterium",
            "dubnium",
            "dysprosium",
            "einsteinium",
            "erbium",
            "europium",
            "fermium",
            "fluorine",
            "francium",
            "gadolinium",
            "gallium",
            "germanium",
            "gold",
            "hafnium",
            "hassium",
            "helium",
            "holmium",
            "hydrogen",
            "indium",
            "iodine",
            "iridium",
            "iron",
            "krypton",
            "lanthanum",
            "lawrencium",
            "lead",
            "lithium",
            "lutetium",
            "magnesium",
            "manganese",
            "meitnerium",
            "mendelevium",
            "mercury",
            "molybdenum",
            "neodymium",
            "neon",
            "neptunium",
            "nickel",
            "niobium",
            "nitrogen",
            "nobelium",
            "osmium",
            "oxygen",
            "palladium",
            "phosphorus",
            "platinum",
            "plutonium",
            "polonium",
            "potassium",
            "praseodymium",
            "promethium",
            "protactinium",
            "radium",
            "radon",
            "rhenium",
            "rhodium",
            "rubidium",
            "ruthenium",
            "rutherfordium",
            "samarium",
            "scandium",
            "seaborgium",
            "selenium",
            "silicon",
            "silver",
            "sodium",
            "strontium",
            "sulfur",
            "tantalum",
            "technetium",
            "tellurium",
            "terbium",
            "thallium",
            "thorium",
            "thulium",
            "tin",
            "titanium",
            "tungsten",
            "uranium",
            "vanadium",
            "xenon",
            "ytterbium",
            "yttrium",
            "zinc",
            "zirconium",
        ] = "actinium",
    ):
        super().__init__(
            **{
                "Category": category,
                "Reds": reds,
                "Greens": greens,
                "Blues": blues,
                "Yellows": yellows,
                "Magentas": magentas,
                "Cyans": cyans,
                "Oranges": oranges,
                "Tints": tints,
                "Grays": grays,
                "Elements": elements,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        category = tree.inputs.menu(
            "Category", "Reds", description="Group of colours, as in PyMOL's Color menu"
        )
        reds = tree.inputs.menu(
            "Reds", "red", description="PyMOL colour, used when Category is Reds"
        )
        greens = tree.inputs.menu(
            "Greens", "green", description="PyMOL colour, used when Category is Greens"
        )
        blues = tree.inputs.menu(
            "Blues", "blue", description="PyMOL colour, used when Category is Blues"
        )
        yellows = tree.inputs.menu(
            "Yellows",
            "yellow",
            description="PyMOL colour, used when Category is Yellows",
        )
        magentas = tree.inputs.menu(
            "Magentas",
            "magenta",
            description="PyMOL colour, used when Category is Magentas",
        )
        cyans = tree.inputs.menu(
            "Cyans", "cyan", description="PyMOL colour, used when Category is Cyans"
        )
        oranges = tree.inputs.menu(
            "Oranges",
            "orange",
            description="PyMOL colour, used when Category is Oranges",
        )
        tints = tree.inputs.menu(
            "Tints", "wheat", description="PyMOL colour, used when Category is Tints"
        )
        grays = tree.inputs.menu(
            "Grays", "white", description="PyMOL colour, used when Category is Grays"
        )
        elements = tree.inputs.menu(
            "Elements",
            "actinium",
            description="PyMOL colour, used when Category is Elements",
        )
        color = tree.outputs.color(
            "Color", (1.0, 1.0, 1.0, 1.0), description="The named colour, in linear RGB"
        )

        menu_switch = g.MenuSwitch.color(
            reds,
            {
                "red": ((1.0, 0.0, 0.0, 1.0), "PyMOL red, #FF0000"),
                "tv_red": ((1.0, 0.033105, 0.033105, 1.0), "PyMOL tv_red, #FF3333"),
                "raspberry": (
                    (0.447988, 0.073239, 0.132868, 1.0),
                    "PyMOL raspberry, #B24C66",
                ),
                "darksalmon": (
                    (0.491905, 0.263273, 0.233022, 1.0),
                    "PyMOL darksalmon, #BA8C85",
                ),
                "salmon": ((1.0, 0.318547, 0.318547, 1.0), "PyMOL salmon, #FF9999"),
                "deepsalmon": (
                    (1.0, 0.214041, 0.214041, 1.0),
                    "PyMOL deepsalmon, #FF8080",
                ),
                "warmpink": (
                    (0.692071, 0.033105, 0.214041, 1.0),
                    "PyMOL warmpink, #D93380",
                ),
                "firebrick": (
                    (0.445146, 0.015325, 0.015325, 1.0),
                    "PyMOL firebrick, #B22121",
                ),
                "ruby": ((0.318547, 0.033105, 0.033105, 1.0), "PyMOL ruby, #993333"),
                "chocolate": (
                    (0.268526, 0.040378, 0.011816, 1.0),
                    "PyMOL chocolate, #8E391C",
                ),
                "brown": ((0.380056, 0.083535, 0.024515, 1.0), "PyMOL brown, #A6522B"),
            },
        )
        menu_switch_1 = g.MenuSwitch.color(
            greens,
            {
                "green": ((0.0, 1.0, 0.0, 1.0), "PyMOL green, #00FF00"),
                "tv_green": ((0.033105, 1.0, 0.033105, 1.0), "PyMOL tv_green, #33FF33"),
                "chartreuse": ((0.214041, 1.0, 0.0, 1.0), "PyMOL chartreuse, #80FF00"),
                "splitpea": ((0.233022, 0.522522, 0.0, 1.0), "PyMOL splitpea, #85BF00"),
                "smudge": (
                    (0.263273, 0.447988, 0.132868, 1.0),
                    "PyMOL smudge, #8CB266",
                ),
                "palegreen": (
                    (0.380056, 0.787412, 0.380056, 1.0),
                    "PyMOL palegreen, #A6E6A6",
                ),
                "limegreen": ((0.0, 1.0, 0.214041, 1.0), "PyMOL limegreen, #00FF80"),
                "lime": ((0.214041, 1.0, 0.214041, 1.0), "PyMOL lime, #80FF80"),
                "limon": ((0.522522, 1.0, 0.050876, 1.0), "PyMOL limon, #BFFF40"),
                "forest": (
                    (0.033105, 0.318547, 0.033105, 1.0),
                    "PyMOL forest, #339933",
                ),
            },
        )
        menu_switch_2 = g.MenuSwitch.color(
            blues,
            {
                "blue": ((0.0, 0.0, 1.0, 1.0), "PyMOL blue, #0000FF"),
                "tv_blue": ((0.073239, 0.073239, 1.0, 1.0), "PyMOL tv_blue, #4C4CFF"),
                "marine": ((0.0, 0.214041, 1.0, 1.0), "PyMOL marine, #0080FF"),
                "slate": ((0.214041, 0.214041, 1.0, 1.0), "PyMOL slate, #8080FF"),
                "lightblue": (
                    (0.522522, 0.522522, 1.0, 1.0),
                    "PyMOL lightblue, #BFBFFF",
                ),
                "skyblue": (
                    (0.033105, 0.214041, 0.603827, 1.0),
                    "PyMOL skyblue, #3380CC",
                ),
                "purpleblue": ((0.214041, 0.0, 1.0, 1.0), "PyMOL purpleblue, #8000FF"),
                "deepblue": (
                    (0.050876, 0.050876, 0.380056, 1.0),
                    "PyMOL deepblue, #4040A6",
                ),
                "density": (
                    (0.010023, 0.010023, 0.318547, 1.0),
                    "PyMOL density, #1A1A99",
                ),
            },
        )
        menu_switch_3 = g.MenuSwitch.color(
            yellows,
            {
                "yellow": ((1.0, 1.0, 0.0, 1.0), "PyMOL yellow, #FFFF00"),
                "tv_yellow": ((1.0, 1.0, 0.033105, 1.0), "PyMOL tv_yellow, #FFFF33"),
                "paleyellow": ((1.0, 1.0, 0.214041, 1.0), "PyMOL paleyellow, #FFFF80"),
                "yelloworange": (
                    (1.0, 0.729347, 0.112805, 1.0),
                    "PyMOL yelloworange, #FFDE5E",
                ),
                "limon": ((0.522522, 1.0, 0.050876, 1.0), "PyMOL limon, #BFFF40"),
                "wheat": ((0.977402, 0.638283, 0.380056, 1.0), "PyMOL wheat, #FCD1A6"),
                "sand": ((0.477, 0.263273, 0.073239, 1.0), "PyMOL sand, #B88C4C"),
            },
        )
        menu_switch_4 = g.MenuSwitch.color(
            magentas,
            {
                "magenta": ((1.0, 0.0, 1.0, 1.0), "PyMOL magenta, #FF00FF"),
                "lightmagenta": (
                    (1.0, 0.033105, 0.603827, 1.0),
                    "PyMOL lightmagenta, #FF33CC",
                ),
                "hotpink": ((1.0, 0.0, 0.214041, 1.0), "PyMOL hotpink, #FF0080"),
                "pink": ((1.0, 0.380056, 0.692071, 1.0), "PyMOL pink, #FFA6D9"),
                "lightpink": (
                    (1.0, 0.522522, 0.729347, 1.0),
                    "PyMOL lightpink, #FFBFDE",
                ),
                "dirtyviolet": (
                    (0.447988, 0.214041, 0.214041, 1.0),
                    "PyMOL dirtyviolet, #B28080",
                ),
                "violet": ((1.0, 0.214041, 1.0, 1.0), "PyMOL violet, #FF80FF"),
                "violetpurple": (
                    (0.263273, 0.050876, 0.318547, 1.0),
                    "PyMOL violetpurple, #8C4099",
                ),
                "purple": ((0.522522, 0.0, 0.522522, 1.0), "PyMOL purple, #BF00BF"),
                "deeppurple": (
                    (0.318547, 0.010023, 0.318547, 1.0),
                    "PyMOL deeppurple, #991A99",
                ),
            },
        )
        menu_switch_5 = g.MenuSwitch.color(
            cyans,
            {
                "cyan": ((0.0, 1.0, 1.0, 1.0), "PyMOL cyan, #00FFFF"),
                "palecyan": ((0.603827, 1.0, 1.0, 1.0), "PyMOL palecyan, #CCFFFF"),
                "aquamarine": ((0.214041, 1.0, 1.0, 1.0), "PyMOL aquamarine, #80FFFF"),
                "greencyan": (
                    (0.050876, 1.0, 0.522522, 1.0),
                    "PyMOL greencyan, #40FFBF",
                ),
                "teal": ((0.0, 0.522522, 0.522522, 1.0), "PyMOL teal, #00BFBF"),
                "deepteal": (
                    (0.010023, 0.318547, 0.318547, 1.0),
                    "PyMOL deepteal, #1A9999",
                ),
                "lightteal": (
                    (0.132868, 0.447988, 0.447988, 1.0),
                    "PyMOL lightteal, #66B2B2",
                ),
            },
        )
        menu_switch_6 = g.MenuSwitch.color(
            oranges,
            {
                "orange": ((1.0, 0.214041, 0.0, 1.0), "PyMOL orange, #FF8000"),
                "tv_orange": (
                    (1.0, 0.263273, 0.019607, 1.0),
                    "PyMOL tv_orange, #FF8C26",
                ),
                "brightorange": (
                    (1.0, 0.447988, 0.033105, 1.0),
                    "PyMOL brightorange, #FFB233",
                ),
                "lightorange": (
                    (1.0, 0.603827, 0.214041, 1.0),
                    "PyMOL lightorange, #FFCC80",
                ),
                "yelloworange": (
                    (1.0, 0.729347, 0.112805, 1.0),
                    "PyMOL yelloworange, #FFDE5E",
                ),
                "olive": ((0.554222, 0.447988, 0.0, 1.0), "PyMOL olive, #C4B200"),
                "deepolive": (
                    (0.318547, 0.318547, 0.010023, 1.0),
                    "PyMOL deepolive, #99991A",
                ),
            },
        )
        menu_switch_7 = g.MenuSwitch.color(
            tints,
            {
                "wheat": ((0.977402, 0.638283, 0.380056, 1.0), "PyMOL wheat, #FCD1A6"),
                "palegreen": (
                    (0.380056, 0.787412, 0.380056, 1.0),
                    "PyMOL palegreen, #A6E6A6",
                ),
                "lightblue": (
                    (0.522522, 0.522522, 1.0, 1.0),
                    "PyMOL lightblue, #BFBFFF",
                ),
                "paleyellow": ((1.0, 1.0, 0.214041, 1.0), "PyMOL paleyellow, #FFFF80"),
                "lightpink": (
                    (1.0, 0.522522, 0.729347, 1.0),
                    "PyMOL lightpink, #FFBFDE",
                ),
                "palecyan": ((0.603827, 1.0, 1.0, 1.0), "PyMOL palecyan, #CCFFFF"),
                "lightorange": (
                    (1.0, 0.603827, 0.214041, 1.0),
                    "PyMOL lightorange, #FFCC80",
                ),
                "bluewhite": (
                    (0.692071, 0.692071, 1.0, 1.0),
                    "PyMOL bluewhite, #D9D9FF",
                ),
            },
        )
        menu_switch_8 = g.MenuSwitch.color(
            grays,
            {
                "white": ((1.0, 1.0, 1.0, 1.0), "PyMOL white, #FFFFFF"),
                "gray90": (
                    (0.805522, 0.805522, 0.805522, 1.0),
                    "PyMOL gray90, #E8E8E8",
                ),
                "gray80": (
                    (0.617615, 0.617615, 0.617615, 1.0),
                    "PyMOL gray80, #CECECE",
                ),
                "gray70": (
                    (0.458124, 0.458124, 0.458124, 1.0),
                    "PyMOL gray70, #B4B4B4",
                ),
                "gray60": (
                    (0.325667, 0.325667, 0.325667, 1.0),
                    "PyMOL gray60, #9B9B9B",
                ),
                "gray50": (
                    (0.218746, 0.218746, 0.218746, 1.0),
                    "PyMOL gray50, #818181",
                ),
                "gray40": (
                    (0.135718, 0.135718, 0.135718, 1.0),
                    "PyMOL gray40, #676767",
                ),
                "gray30": (
                    (0.074748, 0.074748, 0.074748, 1.0),
                    "PyMOL gray30, #4D4D4D",
                ),
                "gray20": (
                    (0.033738, 0.033738, 0.033738, 1.0),
                    "PyMOL gray20, #343434",
                ),
                "gray10": ((0.01018, 0.01018, 0.01018, 1.0), "PyMOL gray10, #1A1A1A"),
                "black": ((0.0, 0.0, 0.0, 1.0), "PyMOL black, #000000"),
            },
        )
        menu_switch_9 = g.MenuSwitch.color(
            elements,
            {
                "actinium": (
                    (0.162029, 0.40724, 0.955973, 1.0),
                    "PyMOL actinium, #70ABFA",
                ),
                "aluminum": (
                    (0.520996, 0.381326, 0.381326, 1.0),
                    "PyMOL aluminum, #BFA6A6",
                ),
                "americium": (
                    (0.088656, 0.107023, 0.887923, 1.0),
                    "PyMOL americium, #545CF2",
                ),
                "antimony": (
                    (0.341914, 0.124772, 0.462077, 1.0),
                    "PyMOL antimony, #9E63B5",
                ),
                "argon": ((0.21586, 0.637597, 0.768151, 1.0), "PyMOL argon, #80D1E3"),
                "arsenic": (
                    (0.508881, 0.21586, 0.768151, 1.0),
                    "PyMOL arsenic, #BD80E3",
                ),
                "astatine": (
                    (0.177888, 0.078187, 0.059511, 1.0),
                    "PyMOL astatine, #754F45",
                ),
                "barium": ((0.0, 0.584078, 0.0, 1.0), "PyMOL barium, #00C900"),
                "berkelium": (
                    (0.254152, 0.078187, 0.768151, 1.0),
                    "PyMOL berkelium, #8A4FE3",
                ),
                "beryllium": ((0.539479, 1.0, 0.0, 1.0), "PyMOL beryllium, #C2FF00"),
                "bismuth": (
                    (0.341914, 0.078187, 0.462077, 1.0),
                    "PyMOL bismuth, #9E4FB5",
                ),
                "bohrium": ((0.745404, 0.0, 0.039546, 1.0), "PyMOL bohrium, #E00038"),
                "boron": ((1.0, 0.462077, 0.462077, 1.0), "PyMOL boron, #FFB5B5"),
                "bromine": (
                    (0.381326, 0.022174, 0.022174, 1.0),
                    "PyMOL bromine, #A62929",
                ),
                "cadmium": ((1.0, 0.693872, 0.274677, 1.0), "PyMOL cadmium, #FFD98F"),
                "calcium": ((0.046665, 1.0, 0.0, 1.0), "PyMOL calcium, #3DFF00"),
                "californium": (
                    (0.3564, 0.036889, 0.658375, 1.0),
                    "PyMOL californium, #A136D4",
                ),
                "carbon": ((0.033105, 1.0, 0.033105, 1.0), "PyMOL carbon, #33FF33"),
                "cerium": ((1.0, 1.0, 0.571125, 1.0), "PyMOL cerium, #FFFFC7"),
                "cesium": (
                    (0.095307, 0.008568, 0.274677, 1.0),
                    "PyMOL cesium, #57178F",
                ),
                "chlorine": (
                    (0.013702, 0.871367, 0.013702, 1.0),
                    "PyMOL chlorine, #1FF01F",
                ),
                "chromium": (
                    (0.254152, 0.318547, 0.571125, 1.0),
                    "PyMOL chromium, #8A99C7",
                ),
                "cobalt": (
                    (0.871367, 0.278894, 0.351533, 1.0),
                    "PyMOL cobalt, #F090A0",
                ),
                "copper": ((0.57758, 0.21586, 0.033105, 1.0), "PyMOL copper, #C88033"),
                "curium": (
                    (0.187821, 0.107023, 0.768151, 1.0),
                    "PyMOL curium, #785CE3",
                ),
                "deuterium": (
                    (0.787412, 0.787412, 0.787412, 1.0),
                    "PyMOL deuterium, #E6E6E6",
                ),
                "dubnium": ((0.637597, 0.0, 0.078187, 1.0), "PyMOL dubnium, #D1004F"),
                "dysprosium": (
                    (0.013702, 1.0, 0.571125, 1.0),
                    "PyMOL dysprosium, #1FFFC7",
                ),
                "einsteinium": (
                    (0.450786, 0.013702, 0.658375, 1.0),
                    "PyMOL einsteinium, #B31FD4",
                ),
                "erbium": ((0.0, 0.791298, 0.177888, 1.0), "PyMOL erbium, #00E675"),
                "europium": ((0.119538, 1.0, 0.571125, 1.0), "PyMOL europium, #61FFC7"),
                "fermium": (
                    (0.450786, 0.013702, 0.491021, 1.0),
                    "PyMOL fermium, #B31FBA",
                ),
                "fluorine": ((0.450786, 1.0, 1.0, 1.0), "PyMOL fluorine, #B3FFFF"),
                "francium": ((0.05448, 0.0, 0.132868, 1.0), "PyMOL francium, #420066"),
                "gadolinium": (
                    (0.059511, 1.0, 0.571125, 1.0),
                    "PyMOL gadolinium, #45FFC7",
                ),
                "gallium": (
                    (0.539479, 0.274677, 0.274677, 1.0),
                    "PyMOL gallium, #C28F8F",
                ),
                "germanium": (
                    (0.132868, 0.274677, 0.274677, 1.0),
                    "PyMOL germanium, #668F8F",
                ),
                "gold": ((1.0, 0.637597, 0.016807, 1.0), "PyMOL gold, #FFD123"),
                "hafnium": ((0.074214, 0.539479, 1.0, 1.0), "PyMOL hafnium, #4DC2FF"),
                "hassium": ((0.791298, 0.0, 0.027321, 1.0), "PyMOL hassium, #E6002E"),
                "helium": ((0.693872, 1.0, 1.0, 1.0), "PyMOL helium, #D9FFFF"),
                "holmium": ((0.0, 1.0, 0.332452, 1.0), "PyMOL holmium, #00FF9C"),
                "hydrogen": (
                    (0.787412, 0.787412, 0.787412, 1.0),
                    "PyMOL hydrogen, #E6E6E6",
                ),
                "indium": (
                    (0.381326, 0.177888, 0.171441, 1.0),
                    "PyMOL indium, #A67573",
                ),
                "iodine": ((0.296138, 0.0, 0.296138, 1.0), "PyMOL iodine, #940094"),
                "iridium": (
                    (0.008568, 0.088656, 0.242281, 1.0),
                    "PyMOL iridium, #175487",
                ),
                "iron": ((0.745404, 0.132868, 0.033105, 1.0), "PyMOL iron, #E06633"),
                "krypton": (
                    (0.107023, 0.47932, 0.637597, 1.0),
                    "PyMOL krypton, #5CB8D1",
                ),
                "lanthanum": (
                    (0.162029, 0.658375, 1.0, 1.0),
                    "PyMOL lanthanum, #70D4FF",
                ),
                "lawrencium": (
                    (0.571125, 0.0, 0.132868, 1.0),
                    "PyMOL lawrencium, #C70066",
                ),
                "lead": ((0.095307, 0.099899, 0.119538, 1.0), "PyMOL lead, #575961"),
                "lithium": ((0.603827, 0.21586, 1.0, 1.0), "PyMOL lithium, #CC80FF"),
                "lutetium": ((0.0, 0.40724, 0.017642, 1.0), "PyMOL lutetium, #00AB24"),
                "magnesium": ((0.254152, 1.0, 0.0, 1.0), "PyMOL magnesium, #8AFF00"),
                "manganese": (
                    (0.332452, 0.194618, 0.571125, 1.0),
                    "PyMOL manganese, #9C7AC7",
                ),
                "meitnerium": (
                    (0.83077, 0.0, 0.019382, 1.0),
                    "PyMOL meitnerium, #EB0026",
                ),
                "mendelevium": (
                    (0.450786, 0.004025, 0.381326, 1.0),
                    "PyMOL mendelevium, #B30DA6",
                ),
                "mercury": (
                    (0.47932, 0.47932, 0.630757, 1.0),
                    "PyMOL mercury, #B8B8D0",
                ),
                "molybdenum": (
                    (0.088656, 0.462077, 0.462077, 1.0),
                    "PyMOL molybdenum, #54B5B5",
                ),
                "neodymium": (
                    (0.571125, 1.0, 0.571125, 1.0),
                    "PyMOL neodymium, #C7FFC7",
                ),
                "neon": ((0.450786, 0.768151, 0.913099, 1.0), "PyMOL neon, #B3E3F5"),
                "neptunium": ((0.0, 0.21586, 1.0, 1.0), "PyMOL neptunium, #0080FF"),
                "nickel": ((0.08022, 0.630757, 0.08022, 1.0), "PyMOL nickel, #50D050"),
                "niobium": (
                    (0.171441, 0.539479, 0.584078, 1.0),
                    "PyMOL niobium, #73C2C9",
                ),
                "nitrogen": ((0.033105, 0.033105, 1.0, 1.0), "PyMOL nitrogen, #3333FF"),
                "nobelium": (
                    (0.508881, 0.004025, 0.242281, 1.0),
                    "PyMOL nobelium, #BD0D87",
                ),
                "osmium": (
                    (0.019382, 0.132868, 0.304987, 1.0),
                    "PyMOL osmium, #266696",
                ),
                "oxygen": ((1.0, 0.073239, 0.073239, 1.0), "PyMOL oxygen, #FF4C4C"),
                "palladium": (
                    (0.0, 0.141263, 0.234551, 1.0),
                    "PyMOL palladium, #006985",
                ),
                "phosphorus": ((1.0, 0.21586, 0.0, 1.0), "PyMOL phosphorus, #FF8000"),
                "platinum": (
                    (0.630757, 0.630757, 0.745404, 1.0),
                    "PyMOL platinum, #D0D0E0",
                ),
                "plutonium": ((0.0, 0.147027, 1.0, 1.0), "PyMOL plutonium, #006BFF"),
                "polonium": ((0.40724, 0.107023, 0.0, 1.0), "PyMOL polonium, #AB5C00"),
                "potassium": (
                    (0.274677, 0.051269, 0.658375, 1.0),
                    "PyMOL potassium, #8F40D4",
                ),
                "praseodymium": (
                    (0.693872, 1.0, 0.571125, 1.0),
                    "PyMOL praseodymium, #D9FFC7",
                ),
                "promethium": (
                    (0.366253, 1.0, 0.571125, 1.0),
                    "PyMOL promethium, #A3FFC7",
                ),
                "protactinium": (
                    (0.0, 0.3564, 1.0, 1.0),
                    "PyMOL protactinium, #00A1FF",
                ),
                "radium": ((0.0, 0.205079, 0.0, 1.0), "PyMOL radium, #007D00"),
                "radon": ((0.05448, 0.223228, 0.304987, 1.0), "PyMOL radon, #428296"),
                "rhenium": (
                    (0.019382, 0.205079, 0.40724, 1.0),
                    "PyMOL rhenium, #267DAB",
                ),
                "rhodium": (
                    (0.003035, 0.205079, 0.262251, 1.0),
                    "PyMOL rhodium, #0A7D8C",
                ),
                "rubidium": (
                    (0.162029, 0.027321, 0.434154, 1.0),
                    "PyMOL rubidium, #702EB0",
                ),
                "ruthenium": (
                    (0.017642, 0.274677, 0.274677, 1.0),
                    "PyMOL ruthenium, #248F8F",
                ),
                "rutherfordium": (
                    (0.603827, 0.0, 0.099899, 1.0),
                    "PyMOL rutherfordium, #CC0059",
                ),
                "samarium": ((0.274677, 1.0, 0.571125, 1.0), "PyMOL samarium, #8FFFC7"),
                "scandium": (
                    (0.791298, 0.791298, 0.791298, 1.0),
                    "PyMOL scandium, #E6E6E6",
                ),
                "seaborgium": (
                    (0.693872, 0.0, 0.059511, 1.0),
                    "PyMOL seaborgium, #D90045",
                ),
                "selenium": ((1.0, 0.3564, 0.0, 1.0), "PyMOL selenium, #FFA100"),
                "silicon": (
                    (0.871367, 0.57758, 0.351533, 1.0),
                    "PyMOL silicon, #F0C8A0",
                ),
                "silver": (
                    (0.527115, 0.527115, 0.527115, 1.0),
                    "PyMOL silver, #C0C0C0",
                ),
                "sodium": ((0.40724, 0.107023, 0.887923, 1.0), "PyMOL sodium, #AB5CF2"),
                "strontium": ((0.0, 1.0, 0.0, 1.0), "PyMOL strontium, #00FF00"),
                "sulfur": (
                    (0.787412, 0.562317, 0.050876, 1.0),
                    "PyMOL sulfur, #E6C640",
                ),
                "tantalum": ((0.074214, 0.381326, 1.0, 1.0), "PyMOL tantalum, #4DA6FF"),
                "technetium": (
                    (0.043735, 0.341914, 0.341914, 1.0),
                    "PyMOL technetium, #3B9E9E",
                ),
                "tellurium": (
                    (0.658375, 0.194618, 0.0, 1.0),
                    "PyMOL tellurium, #D47A00",
                ),
                "terbium": ((0.029557, 1.0, 0.571125, 1.0), "PyMOL terbium, #30FFC7"),
                "thallium": (
                    (0.381326, 0.088656, 0.074214, 1.0),
                    "PyMOL thallium, #A6544D",
                ),
                "thorium": ((0.0, 0.491021, 1.0, 1.0), "PyMOL thorium, #00BAFF"),
                "thulium": ((0.0, 0.658375, 0.084376, 1.0), "PyMOL thulium, #00D452"),
                "tin": ((0.132868, 0.21586, 0.21586, 1.0), "PyMOL tin, #668080"),
                "titanium": (
                    (0.520996, 0.539479, 0.571125, 1.0),
                    "PyMOL titanium, #BFC2C7",
                ),
                "tungsten": (
                    (0.015209, 0.296138, 0.672443, 1.0),
                    "PyMOL tungsten, #2194D6",
                ),
                "uranium": ((0.0, 0.274677, 1.0, 1.0), "PyMOL uranium, #008FFF"),
                "vanadium": (
                    (0.381326, 0.381326, 0.40724, 1.0),
                    "PyMOL vanadium, #A6A6AB",
                ),
                "xenon": ((0.05448, 0.341914, 0.434154, 1.0), "PyMOL xenon, #429EB0"),
                "ytterbium": (
                    (0.0, 0.520996, 0.039546, 1.0),
                    "PyMOL ytterbium, #00BF38",
                ),
                "yttrium": ((0.296138, 1.0, 1.0, 1.0), "PyMOL yttrium, #94FFFF"),
                "zinc": ((0.205079, 0.21586, 0.434154, 1.0), "PyMOL zinc, #7D80B0"),
                "zirconium": (
                    (0.296138, 0.745404, 0.745404, 1.0),
                    "PyMOL zirconium, #94E0E0",
                ),
            },
        )
        (
            g.MenuSwitch.color(
                category,
                {
                    "Reds": menu_switch.o.output,
                    "Greens": menu_switch_1.o.output,
                    "Blues": menu_switch_2.o.output,
                    "Yellows": menu_switch_3.o.output,
                    "Magentas": menu_switch_4.o.output,
                    "Cyans": menu_switch_5.o.output,
                    "Oranges": menu_switch_6.o.output,
                    "Tints": menu_switch_7.o.output,
                    "Grays": menu_switch_8.o.output,
                    "Elements": menu_switch_9.o.output,
                },
            )
            >> color
        )


ASSET = ColorPyMOL

ASSET_METADATA = {
    "description": "A named colour from PyMOL, including its element colours",
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
