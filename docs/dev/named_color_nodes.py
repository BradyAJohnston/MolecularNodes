"""
Generate the "Color PyMOL", "Color VMD" and "Color CSS" nodes from the named
colours that `molecularnodes.color.parse_color` reads.

    uv run docs/dev/named_color_nodes.py
    uv run -m nodebpy.assets build && uv run -m nodebpy.assets dump

PyMOL and VMD colours come from `molecularnodes/_color_names.py`, CSS colours
(MolViewSpec, ChimeraX) from matplotlib, so the nodes and `parse_color` agree.
Each node is a Menu Switch of constant colours, converted from sRGB to linear.
The PyMOL node groups its colours as PyMOL's own Color menu does
(`all_colors_list` in modules/pymol/menu.py), with the element colours as one
more category.
"""

import importlib.util
from pathlib import Path
import numpy as np
from matplotlib.colors import CSS4_COLORS, to_rgb

ROOT = Path(__file__).parents[2]
GEOMETRY = ROOT / "molecularnodes/nodes/geometry"
CATALOG_ID = "d3f975df-8408-4972-a669-8187a57e01d0"  # Color

# loaded from the file so that generating doesn't need the add-on or bpy
_spec = importlib.util.spec_from_file_location(
    "_color_names", ROOT / "molecularnodes/_color_names.py"
)
_names = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_names)
PYMOL, VMD = _names.PYMOL_COLORS, _names.VMD_COLORS

# PyMOL's Color menu, from all_colors_list in modules/pymol/menu.py
PYMOL_MENU = {
    "Reds": [
        "red", "tv_red", "raspberry", "darksalmon", "salmon", "deepsalmon",
        "warmpink", "firebrick", "ruby", "chocolate", "brown",
    ],
    "Greens": [
        "green", "tv_green", "chartreuse", "splitpea", "smudge", "palegreen",
        "limegreen", "lime", "limon", "forest",
    ],
    "Blues": [
        "blue", "tv_blue", "marine", "slate", "lightblue", "skyblue",
        "purpleblue", "deepblue", "density",
    ],
    "Yellows": [
        "yellow", "tv_yellow", "paleyellow", "yelloworange", "limon", "wheat",
        "sand",
    ],
    "Magentas": [
        "magenta", "lightmagenta", "hotpink", "pink", "lightpink", "dirtyviolet",
        "violet", "violetpurple", "purple", "deeppurple",
    ],
    "Cyans": [
        "cyan", "palecyan", "aquamarine", "greencyan", "teal", "deepteal",
        "lightteal",
    ],
    "Oranges": [
        "orange", "tv_orange", "brightorange", "lightorange", "yelloworange",
        "olive", "deepolive",
    ],
    "Tints": [
        "wheat", "palegreen", "lightblue", "paleyellow", "lightpink", "palecyan",
        "lightorange", "bluewhite",
    ],
    "Grays": [
        "white", "gray90", "gray80", "gray70", "gray60", "gray50", "gray40",
        "gray30", "gray20", "gray10", "black",
    ],
}  # fmt: skip
_order = list(PYMOL)
PYMOL_MENU["Elements"] = sorted(
    ["hydrogen", "deuterium", "carbon", "nitrogen", "oxygen", "sulfur"]
    + _order[_order.index("helium") : _order.index("meitnerium") + 1]
)


def to_linear(c):
    c = np.asarray(c, dtype=float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def to_hex(rgb) -> str:
    return "#" + "".join(f"{round(v * 255):02X}" for v in rgb)


def items(colors: dict[str, tuple], label: str) -> str:
    """The Menu Switch items dict: name -> (linear RGBA, description)."""
    lines = []
    for name, rgb in colors.items():
        r, g, b = (round(float(v), 6) for v in to_linear(rgb))
        description = f"{label} {name}, {to_hex(rgb)}"
        lines.append(f'"{name}": (({r}, {g}, {b}, 1.0), "{description}"),')
    return "{" + "\n".join(lines) + "}"


HEADER = '''from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import AssetGeometryGroup, PackageLibrary


class {cls}(AssetGeometryGroup):
    """{name}"""

    _name = "{name}"
    _asset_name = "{name}"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
'''

FOOTER = """

ASSET = {cls}

ASSET_METADATA = {{
    "catalog_id": "{catalog}",
    "description": "{description}",
}}
"""


def single_menu(cls, name, colors, default, label, menu_description, description):
    return (
        HEADER.format(cls=cls, name=name)
        + f'''        menu = tree.inputs.menu("Color", "{default}", description="{menu_description}")
        color = tree.outputs.color("Color", description="The named colour, in linear RGB")
        g.MenuSwitch.color(menu, {items(colors, label)}) >> color
'''
        + FOOTER.format(cls=cls, catalog=CATALOG_ID, description=description)
    )


def pymol() -> str:
    body = [
        '        category = tree.inputs.menu("Category", "Reds", '
        'description="Group of colours, as in PyMOL\'s Color menu")'
    ]
    for group, names in PYMOL_MENU.items():
        var = group.lower()
        body.append(
            f'        {var} = tree.inputs.menu("{group}", "{names[0]}", '
            f'description="PyMOL colour, used when Category is {group}")'
        )
    body.append(
        '        color = tree.outputs.color("Color", '
        'description="The named colour, in linear RGB")'
    )
    switches = []
    for group, names in PYMOL_MENU.items():
        var = group.lower()
        body.append(
            f"        {var}_color = g.MenuSwitch.color({var}, "
            f"{items({n: PYMOL[n] for n in names}, 'PyMOL')})"
        )
        switches.append(f'"{group}": {var}_color.o.output')
    body.append(
        f"        g.MenuSwitch.color(category, {{{', '.join(switches)}}}) >> color"
    )
    return (
        HEADER.format(cls="ColorPyMOL", name="Color PyMOL")
        + "\n".join(body)
        + "\n"
        + FOOTER.format(
            cls="ColorPyMOL",
            catalog=CATALOG_ID,
            description="A named colour from PyMOL, including its element colours",
        )
    )


def vmd() -> str:
    return single_menu(
        "ColorVMD",
        "Color VMD",
        VMD,
        "blue",
        "VMD",
        "VMD colour, in colour ID order (0 to 32)",
        "A named colour from VMD",
    )


def css() -> str:
    colors = {name: to_rgb(value) for name, value in sorted(CSS4_COLORS.items())}
    return single_menu(
        "ColorCSS",
        "Color CSS",
        colors,
        "steelblue",
        "CSS",
        "CSS colour name, as used by MolViewSpec and ChimeraX",
        "A named CSS colour, as used by MolViewSpec and ChimeraX",
    )


if __name__ == "__main__":
    for filename, source in (
        ("color_pymol.py", pymol()),
        ("color_vmd.py", vmd()),
        ("color_css.py", css()),
    ):
        (GEOMETRY / filename).write_text(source)
        print(f"wrote {GEOMETRY / filename}")
