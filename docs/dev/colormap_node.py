"""
Generate the "Color matplotlib" node from matplotlib's colormaps.

    uv run docs/dev/colormap_node.py
    uv run -m nodebpy.assets build && uv run -m nodebpy.assets dump

Writes ``molecularnodes/nodes/geometry/color_matplotlib.py`` with the ramp data
inline; build + dump turns it into the literal tree source that is committed.
matplotlib is needed only here, never at runtime. Re-run when matplotlib adds
colormaps: the script fails if a map is not assigned to a category below.

Each map becomes one Color Ramp of at most 32 stops (Blender's limit, beyond
which extra stops are silently dropped). Stops are chosen greedily where the
linear interpolation is furthest from matplotlib's 256 samples, working in
linear RGB because that is what the ramp stores and interpolates; the error
is measured back in 8-bit sRGB. Small ListedColormaps (the qualitative maps)
are copied exactly with constant interpolation.
"""

from pathlib import Path
import matplotlib
import numpy as np
from matplotlib.colors import ListedColormap

MAX_STOPS = 32
TOLERANCE = 0.5  # stop adding stops once every sample is within this, in 1/255 sRGB
LOSSY = 2.0  # errors above this are called out in the node's frame comment

# matplotlib's own grouping, from its colormap reference gallery. Menu name,
# full name for the tooltip, maps.
CATEGORIES = {
    "Uniform": (
        "Perceptually uniform sequential",
        ["viridis", "plasma", "inferno", "magma", "cividis"],
    ),
    "Sequential": (
        "Sequential",
        [
            "Greys", "Purples", "Blues", "Greens", "Oranges", "Reds", "YlOrBr",
            "YlOrRd", "OrRd", "PuRd", "RdPu", "BuPu", "GnBu", "PuBu", "YlGnBu",
            "PuBuGn", "BuGn", "YlGn",
        ],
    ),
    "Sequential 2": (
        "Sequential (2)",
        [
            "binary", "gray", "bone", "pink", "spring", "summer", "autumn",
            "winter", "cool", "Wistia", "hot", "afmhot", "gist_heat", "copper",
        ],
    ),
    "Diverging": (
        "Diverging",
        [
            "PiYG", "PRGn", "BrBG", "PuOr", "RdGy", "RdBu", "RdYlBu", "RdYlGn",
            "Spectral", "coolwarm", "bwr", "seismic", "berlin", "managua", "vanimo",
        ],
    ),
    "Cyclic": ("Cyclic", ["twilight", "twilight_shifted", "hsv"]),
    "Qualitative": (
        "Qualitative",
        [
            "Pastel1", "Pastel2", "Paired", "Accent", "Dark2", "Set1", "Set2",
            "Set3", "tab10", "tab20", "tab20b", "tab20c", "okabe_ito",
        ],
    ),
    "Miscellaneous": (
        "Miscellaneous",
        [
            "ocean", "gist_earth", "terrain", "gist_stern", "gnuplot", "gnuplot2",
            "CMRmap", "cubehelix", "brg", "gist_rainbow", "rainbow", "jet",
            "turbo", "nipy_spectral", "gist_ncar",
        ],
    ),
}  # fmt: skip

# high-frequency stripes that no 32-stop ramp can represent
EXCLUDED = {"flag", "prism"}

OUT = Path(__file__).parents[2] / "molecularnodes/nodes/geometry/color_matplotlib.py"


def to_linear(c):
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def to_srgb(c):
    c = np.clip(c, 0.0, None)
    return np.where(c <= 0.0031308, c * 12.92, 1.055 * c ** (1 / 2.4) - 0.055)


def samples(name: str) -> np.ndarray:
    return matplotlib.colormaps[name](np.linspace(0.0, 1.0, 256))[:, :3]


def fit(name: str):
    """(stops, interpolation, max error in 1/255 sRGB) for one colormap."""
    cmap = matplotlib.colormaps[name]
    if isinstance(cmap, ListedColormap) and cmap.N <= MAX_STOPS:
        # matplotlib maps x to colour floor(x * N); a constant ramp with stop i
        # at i / N does the same
        colors = to_linear(np.asarray(cmap.colors, dtype=float)[:, :3])
        stops = [(i / cmap.N, c) for i, c in enumerate(colors)]
        return stops, "CONSTANT", 0.0
    srgb = samples(name)
    lin = to_linear(srgb)
    x = np.linspace(0.0, 1.0, len(lin))
    chosen = [0, len(lin) - 1]
    while True:
        s = sorted(chosen)
        approx = np.stack([np.interp(x, x[s], lin[s, k]) for k in range(3)], axis=1)
        error = np.abs(to_srgb(approx) - srgb).max(axis=1) * 255
        if len(s) >= MAX_STOPS or error.max() < TOLERANCE:
            return [(x[i], lin[i]) for i in s], "LINEAR", float(error.max())
        chosen.append(int(error.argmax()))


def check_coverage() -> None:
    """Every matplotlib map is categorised, excluded, or an alias of one that is."""
    listed = {m for _, maps in CATEGORIES.values() for m in maps}
    for name in sorted(matplotlib.colormaps):
        if name.endswith("_r") or name in listed or name in EXCLUDED:
            continue
        twin = next(
            (m for m in listed if np.array_equal(samples(m), samples(name))), None
        )
        if twin is None:
            raise SystemExit(f"{name!r} is not in CATEGORIES; add it or exclude it.")


def comment(name: str, category: str, n: int, interpolation: str, error: float) -> str:
    if interpolation == "CONSTANT":
        text = (
            f"{name}: {category}. {n} colours copied exactly, constant interpolation."
        )
    else:
        text = (
            f"{name}: {category}. {n} stops fitted to matplotlib "
            f"{matplotlib.__version__} in linear RGB, max error {error:.1f}/255 in sRGB."
        )
        if error > LOSSY:
            text += " Detail finer than the stops can hold is lost."
    return text


def ramp_items(stops) -> str:
    return ", ".join(
        f"({round(float(p), 6)}, ({', '.join(str(round(float(v), 6)) for v in c)}, 1.0))"
        for p, c in stops
    )


def generate() -> str:
    lines = []
    for menu, (category, maps) in CATEGORIES.items():
        var = menu.lower().replace(" ", "_")
        lines.append(f"        {var}_items = {{}}")
        for name in maps:
            stops, interpolation, error = fit(name)
            text = comment(name, category, len(stops), interpolation, error)
            lines += [
                f'        with g.Frame("{name}"):',
                f"            g.String(string={text!r})",
                f"            {var}_items[{name!r}] = (",
                f"                g.ColorRamp(value, items=[{ramp_items(stops)}], "
                f'color_interpolation="{interpolation}").o.color,',
                f"                {category + ' colormap from matplotlib'!r},",
                "            )",
            ]
        lines.append(f"        {var}_switch = g.MenuSwitch.color({var}, {var}_items)")
    switches = ", ".join(
        f'"{menu}": ({menu.lower().replace(" ", "_")}_switch.o.output, {category!r})'
        for menu, (category, _) in CATEGORIES.items()
    )
    lines.append(f"        g.MenuSwitch.color(category, {{{switches}}}) >> color")
    lines.append('        category.default_value = "Uniform"')
    for menu, (_, maps) in CATEGORIES.items():
        lines.append(
            f'        {menu.lower().replace(" ", "_")}.default_value = "{maps[0]}"'
        )

    menus = "\n".join(
        f"        {menu.lower().replace(' ', '_')} = tree.inputs.menu("
        f'"{menu}", description="{category} colormap, used when Category is {menu}")'
        for menu, (category, _) in CATEGORIES.items()
    )
    return f'''from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import AssetGeometryGroup, PackageLibrary


class ColorMatplotlib(AssetGeometryGroup):
    """Color matplotlib"""

    _name = "Color matplotlib"
    _asset_name = "Color matplotlib"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"
    _tree_properties = {{"node_tool_idname": "geometry.color_matplotlib"}}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        value = tree.inputs.float(
            "Value",
            0.5,
            description="Position along the colormap, 0 to 1. Values outside are clamped",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        reverse = tree.inputs.boolean(
            "Reverse", False, description="Run the colormap from 1 to 0, like matplotlib's _r maps"
        )
        category = tree.inputs.menu(
            "Category", description="Group of colormaps, as in matplotlib's colormap reference"
        )
{menus}
        color = tree.outputs.color(
            "Color", description="Colour of the selected colormap at Value"
        )

        value = g.Switch.float(reverse, value, 1.0 - value)
{chr(10).join(lines)}


ASSET = ColorMatplotlib

ASSET_METADATA = {{
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
    "description": "Colour a value from 0 to 1 with any of matplotlib's colormaps",
}}
'''


if __name__ == "__main__":
    check_coverage()
    OUT.write_text(generate())
    print(f"wrote {OUT}")
