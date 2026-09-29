# Node-group asset "Color matplotlib" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ColorSocket,
    FloatSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputMenu


class ColorMatplotlib(AssetGeometryGroup):
    """
    Color matplotlib

    Parameters
    ----------
    value : InputFloat
        Position along the colormap, 0 to 1. Values outside are clamped
    reverse : InputBoolean
        Run the colormap from 1 to 0, like matplotlib's _r maps
    category : InputMenu | Literal["Uniform", "Sequential", "Sequential 2", "Diverging", "Cyclic", "Qualitative", "Miscellaneous"]
        Group of colormaps, as in matplotlib's colormap reference
    uniform : InputMenu | Literal["viridis", "plasma", "inferno", "magma", "cividis"]
        Perceptually uniform sequential colormap, used when Category is Uniform
    sequential : InputMenu | Literal["Greys", "Purples", "Blues", "Greens", "Oranges", "Reds", "YlOrBr", "YlOrRd", "OrRd", "PuRd", "RdPu", "BuPu", "GnBu", "PuBu", "YlGnBu", "PuBuGn", "BuGn", "YlGn"]
        Sequential colormap, used when Category is Sequential
    sequential_2 : InputMenu | Literal["binary", "gray", "bone", "pink", "spring", "summer", "autumn", "winter", "cool", "Wistia", "hot", "afmhot", "gist_heat", "copper"]
        Sequential (2) colormap, used when Category is Sequential 2
    diverging : InputMenu | Literal["PiYG", "PRGn", "BrBG", "PuOr", "RdGy", "RdBu", "RdYlBu", "RdYlGn", "Spectral", "coolwarm", "bwr", "seismic", "berlin", "managua", "vanimo"]
        Diverging colormap, used when Category is Diverging
    cyclic : InputMenu | Literal["twilight", "twilight_shifted", "hsv"]
        Cyclic colormap, used when Category is Cyclic
    qualitative : InputMenu | Literal["Pastel1", "Pastel2", "Paired", "Accent", "Dark2", "Set1", "Set2", "Set3", "tab10", "tab20", "tab20b", "tab20c", "okabe_ito"]
        Qualitative colormap, used when Category is Qualitative
    miscellaneous : InputMenu | Literal["ocean", "gist_earth", "terrain", "gist_stern", "gnuplot", "gnuplot2", "CMRmap", "cubehelix", "brg", "gist_rainbow", "rainbow", "jet", "turbo", "nipy_spectral", "gist_ncar"]
        Miscellaneous colormap, used when Category is Miscellaneous

    Inputs
    ------
    i.value : FloatSocket
        Position along the colormap, 0 to 1. Values outside are clamped
    i.reverse : BooleanSocket
        Run the colormap from 1 to 0, like matplotlib's _r maps
    i.category : MenuSocket
        Group of colormaps, as in matplotlib's colormap reference
    i.uniform : MenuSocket
        Perceptually uniform sequential colormap, used when Category is Uniform
    i.sequential : MenuSocket
        Sequential colormap, used when Category is Sequential
    i.sequential_2 : MenuSocket
        Sequential (2) colormap, used when Category is Sequential 2
    i.diverging : MenuSocket
        Diverging colormap, used when Category is Diverging
    i.cyclic : MenuSocket
        Cyclic colormap, used when Category is Cyclic
    i.qualitative : MenuSocket
        Qualitative colormap, used when Category is Qualitative
    i.miscellaneous : MenuSocket
        Miscellaneous colormap, used when Category is Miscellaneous

    Outputs
    -------
    o.color : ColorSocket
        Colour of the selected colormap at Value
    """

    _name = "Color matplotlib"
    _asset_name = "Color matplotlib"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_matplotlib"}

    class _Inputs(SocketAccessor):
        value: FloatSocket
        """Position along the colormap, 0 to 1. Values outside are clamped"""
        reverse: BooleanSocket
        """Run the colormap from 1 to 0, like matplotlib's _r maps"""
        category: MenuSocket
        """Group of colormaps, as in matplotlib's colormap reference"""
        uniform: MenuSocket
        """Perceptually uniform sequential colormap, used when Category is Uniform"""
        sequential: MenuSocket
        """Sequential colormap, used when Category is Sequential"""
        sequential_2: MenuSocket
        """Sequential (2) colormap, used when Category is Sequential 2"""
        diverging: MenuSocket
        """Diverging colormap, used when Category is Diverging"""
        cyclic: MenuSocket
        """Cyclic colormap, used when Category is Cyclic"""
        qualitative: MenuSocket
        """Qualitative colormap, used when Category is Qualitative"""
        miscellaneous: MenuSocket
        """Miscellaneous colormap, used when Category is Miscellaneous"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """Colour of the selected colormap at Value"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputFloat = 0.5,
        reverse: InputBoolean = False,
        category: InputMenu
        | Literal[
            "Uniform",
            "Sequential",
            "Sequential 2",
            "Diverging",
            "Cyclic",
            "Qualitative",
            "Miscellaneous",
        ] = "Uniform",
        uniform: InputMenu
        | Literal["viridis", "plasma", "inferno", "magma", "cividis"] = "viridis",
        sequential: InputMenu
        | Literal[
            "Greys",
            "Purples",
            "Blues",
            "Greens",
            "Oranges",
            "Reds",
            "YlOrBr",
            "YlOrRd",
            "OrRd",
            "PuRd",
            "RdPu",
            "BuPu",
            "GnBu",
            "PuBu",
            "YlGnBu",
            "PuBuGn",
            "BuGn",
            "YlGn",
        ] = "Greys",
        sequential_2: InputMenu
        | Literal[
            "binary",
            "gray",
            "bone",
            "pink",
            "spring",
            "summer",
            "autumn",
            "winter",
            "cool",
            "Wistia",
            "hot",
            "afmhot",
            "gist_heat",
            "copper",
        ] = "binary",
        diverging: InputMenu
        | Literal[
            "PiYG",
            "PRGn",
            "BrBG",
            "PuOr",
            "RdGy",
            "RdBu",
            "RdYlBu",
            "RdYlGn",
            "Spectral",
            "coolwarm",
            "bwr",
            "seismic",
            "berlin",
            "managua",
            "vanimo",
        ] = "PiYG",
        cyclic: InputMenu | Literal["twilight", "twilight_shifted", "hsv"] = "twilight",
        qualitative: InputMenu
        | Literal[
            "Pastel1",
            "Pastel2",
            "Paired",
            "Accent",
            "Dark2",
            "Set1",
            "Set2",
            "Set3",
            "tab10",
            "tab20",
            "tab20b",
            "tab20c",
            "okabe_ito",
        ] = "Pastel1",
        miscellaneous: InputMenu
        | Literal[
            "ocean",
            "gist_earth",
            "terrain",
            "gist_stern",
            "gnuplot",
            "gnuplot2",
            "CMRmap",
            "cubehelix",
            "brg",
            "gist_rainbow",
            "rainbow",
            "jet",
            "turbo",
            "nipy_spectral",
            "gist_ncar",
        ] = "ocean",
    ):
        super().__init__(
            **{
                "Value": value,
                "Reverse": reverse,
                "Category": category,
                "Uniform": uniform,
                "Sequential": sequential,
                "Sequential 2": sequential_2,
                "Diverging": diverging,
                "Cyclic": cyclic,
                "Qualitative": qualitative,
                "Miscellaneous": miscellaneous,
            }
        )

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
            "Reverse",
            False,
            description="Run the colormap from 1 to 0, like matplotlib's _r maps",
        )
        category = tree.inputs.menu(
            "Category",
            description="Group of colormaps, as in matplotlib's colormap reference",
        )
        uniform = tree.inputs.menu(
            "Uniform",
            description="Perceptually uniform sequential colormap, used when Category is Uniform",
        )
        sequential = tree.inputs.menu(
            "Sequential",
            description="Sequential colormap, used when Category is Sequential",
        )
        sequential_2 = tree.inputs.menu(
            "Sequential 2",
            description="Sequential (2) colormap, used when Category is Sequential 2",
        )
        diverging = tree.inputs.menu(
            "Diverging",
            description="Diverging colormap, used when Category is Diverging",
        )
        cyclic = tree.inputs.menu(
            "Cyclic", description="Cyclic colormap, used when Category is Cyclic"
        )
        qualitative = tree.inputs.menu(
            "Qualitative",
            description="Qualitative colormap, used when Category is Qualitative",
        )
        miscellaneous = tree.inputs.menu(
            "Miscellaneous",
            description="Miscellaneous colormap, used when Category is Miscellaneous",
        )
        color = tree.outputs.color(
            "Color",
            (1.0, 1.0, 1.0, 1.0),
            description="Colour of the selected colormap at Value",
        )

        switch = reverse.switch.float(value, 1.0 - value)
        with g.Frame("viridis"):
            _string = g.String(
                string="viridis: Perceptually uniform sequential. 13 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.056964, 0.0, 0.084614, 1.0)),
                    (0.033127, (0.063545, 0.003412, 0.116646, 1.0)),
                    (0.091617, (0.066458, 0.013306, 0.169059, 1.0)),
                    (0.178541, (0.059207, 0.041897, 0.242128, 1.0)),
                    (0.295513, (0.033152, 0.106746, 0.269137, 1.0)),
                    (0.439554, (0.019165, 0.21278, 0.277688, 1.0)),
                    (0.583146, (0.007237, 0.371209, 0.2478, 1.0)),
                    (0.693999, (0.037512, 0.512948, 0.170768, 1.0)),
                    (0.786693, (0.149213, 0.628924, 0.090451, 1.0)),
                    (0.864936, (0.35698, 0.713773, 0.032054, 1.0)),
                    (0.928568, (0.624903, 0.754318, 0.007932, 1.0)),
                    (0.975252, (0.861662, 0.786243, 0.009231, 1.0)),
                    (1.0, (1.0, 0.800834, 0.019509, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("plasma"):
            _string_1 = g.String(
                string="plasma: Perceptually uniform sequential. 12 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_1 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.002266, 0.23066, 1.0)),
                    (0.051716, (0.025273, 0.001613, 0.306928, 1.0)),
                    (0.110725, (0.054617, 0.000757, 0.340053, 1.0)),
                    (0.197836, (0.134493, 0.0, 0.407742, 1.0)),
                    (0.306829, (0.277858, 0.000498, 0.384232, 1.0)),
                    (0.44409, (0.510609, 0.029375, 0.226157, 1.0)),
                    (0.593944, (0.747063, 0.111242, 0.118938, 1.0)),
                    (0.72879, (0.939292, 0.255711, 0.054775, 1.0)),
                    (0.836552, (0.998643, 0.448827, 0.025601, 1.0)),
                    (0.920225, (0.973129, 0.671966, 0.013724, 1.0)),
                    (0.975784, (0.894391, 0.865177, 0.024056, 1.0)),
                    (1.0, (0.865465, 0.954195, 0.014721, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("inferno"):
            _string_2 = g.String(
                string="inferno: Perceptually uniform sequential. 14 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_2 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.00004, 0.000213, 1.0)),
                    (0.044799, (0.001378, 0.000747, 0.006413, 1.0)),
                    (0.10752, (0.00649, 0.005094, 0.038209, 1.0)),
                    (0.190219, (0.036722, 0.001158, 0.149166, 1.0)),
                    (0.308832, (0.140889, 0.00799, 0.16775, 1.0)),
                    (0.44726, (0.358245, 0.023378, 0.120673, 1.0)),
                    (0.589525, (0.721253, 0.058802, 0.040077, 1.0)),
                    (0.709145, (0.930775, 0.183965, 0.003106, 1.0)),
                    (0.797535, (0.987112, 0.365903, 0.0, 1.0)),
                    (0.85714, (0.961455, 0.530481, 0.013629, 1.0)),
                    (0.90899, (0.913782, 0.707313, 0.063495, 1.0)),
                    (0.946232, (0.876533, 0.847098, 0.149092, 1.0)),
                    (0.975791, (0.888149, 0.942092, 0.260341, 1.0)),
                    (1.0, (0.987122, 1.0, 0.391105, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("magma"):
            _string_3 = g.String(
                string="magma: Perceptually uniform sequential. 14 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.7/255 in sRGB."
            )
            color_ramp_3 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.000033, 1.0)),
                    (0.046137, (0.00124, 0.000962, 0.006846, 1.0)),
                    (0.102791, (0.006375, 0.005145, 0.033288, 1.0)),
                    (0.159361, (0.019424, 0.006345, 0.102104, 1.0)),
                    (0.224753, (0.054897, 0.003502, 0.199433, 1.0)),
                    (0.313173, (0.133709, 0.010451, 0.219376, 1.0)),
                    (0.419181, (0.276656, 0.024971, 0.225683, 1.0)),
                    (0.5411, (0.554481, 0.041017, 0.18117, 1.0)),
                    (0.660554, (0.914511, 0.089372, 0.09224, 1.0)),
                    (0.76244, (0.979965, 0.261004, 0.117131, 1.0)),
                    (0.846355, (1.0, 0.454336, 0.192139, 1.0)),
                    (0.917321, (0.987665, 0.680951, 0.316856, 1.0)),
                    (0.96712, (0.976634, 0.847427, 0.425592, 1.0)),
                    (1.0, (0.970055, 1.0, 0.535729, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("cividis"):
            _string_4 = g.String(
                string="cividis: Perceptually uniform sequential. 13 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.7/255 in sRGB."
            )
            color_ramp_4 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.016467, 0.075244, 1.0)),
                    (0.030621, (0.0, 0.021184, 0.103398, 1.0)),
                    (0.079477, (0.0, 0.030102, 0.164417, 1.0)),
                    (0.090894, (0.000187, 0.031563, 0.163847, 1.0)),
                    (0.18086, (0.028134, 0.054017, 0.151467, 1.0)),
                    (0.291146, (0.073707, 0.091234, 0.149138, 1.0)),
                    (0.402322, (0.133844, 0.141348, 0.163683, 1.0)),
                    (0.519735, (0.216215, 0.210928, 0.191532, 1.0)),
                    (0.644381, (0.346179, 0.308209, 0.180823, 1.0)),
                    (0.769231, (0.524198, 0.437177, 0.145954, 1.0)),
                    (0.884914, (0.740223, 0.593282, 0.094374, 1.0)),
                    (0.987296, (0.985118, 0.768256, 0.033705, 1.0)),
                    (1.0, (0.992194, 0.807711, 0.038779, 1.0)),
                ),
            )
        with g.Frame("Greys"):
            _string_5 = g.String(
                string="Greys: Sequential. 10 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_5 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (1.0, 1.0, 1.0, 1.0)),
                    (0.04759, (0.956198, 0.956198, 0.956198, 1.0)),
                    (0.174892, (0.837065, 0.837065, 0.837065, 1.0)),
                    (0.395811, (0.450615, 0.450615, 0.450615, 1.0)),
                    (0.569453, (0.202972, 0.202972, 0.202972, 1.0)),
                    (0.743593, (0.078739, 0.078739, 0.078739, 1.0)),
                    (0.851876, (0.019945, 0.019945, 0.019945, 1.0)),
                    (0.93213, (0.005848, 0.005848, 0.005848, 1.0)),
                    (0.983865, (0.000984, 0.000984, 0.000984, 1.0)),
                    (1.0, (0.0, 0.0, 0.0, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("Purples"):
            _string_6 = g.String(
                string="Purples: Sequential. 11 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_6 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.993746, 0.988751, 1.0, 1.0)),
                    (0.162947, (0.827364, 0.802343, 0.877711, 1.0)),
                    (0.312004, (0.601548, 0.625307, 0.798818, 1.0)),
                    (0.471388, (0.360511, 0.328124, 0.585732, 1.0)),
                    (0.602077, (0.228701, 0.232305, 0.515827, 1.0)),
                    (0.709373, (0.158486, 0.10378, 0.399843, 1.0)),
                    (0.793046, (0.125541, 0.053332, 0.328024, 1.0)),
                    (0.869115, (0.086838, 0.018918, 0.278668, 1.0)),
                    (0.933115, (0.069933, 0.006481, 0.238384, 1.0)),
                    (0.98359, (0.051817, 0.001008, 0.212146, 1.0)),
                    (1.0, (0.049954, 0.0, 0.203507, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("Blues"):
            _string_7 = g.String(
                string="Blues: Sequential. 14 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.7/255 in sRGB."
            )
            color_ramp_7 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.928189, 0.962965, 0.999953, 1.0)),
                    (0.122139, (0.733011, 0.834867, 0.930923, 1.0)),
                    (0.246143, (0.567173, 0.709106, 0.866029, 1.0)),
                    (0.356595, (0.369525, 0.610514, 0.767134, 1.0)),
                    (0.447362, (0.216626, 0.492207, 0.704492, 1.0)),
                    (0.500088, (0.145711, 0.421336, 0.672954, 1.0)),
                    (0.566307, (0.090451, 0.347824, 0.613344, 1.0)),
                    (0.624649, (0.054041, 0.286997, 0.565128, 1.0)),
                    (0.685776, (0.031345, 0.222462, 0.512966, 1.0)),
                    (0.747419, (0.01532, 0.166484, 0.464197, 1.0)),
                    (0.81362, (0.006901, 0.118218, 0.392864, 1.0)),
                    (0.87319, (0.002389, 0.082508, 0.333363, 1.0)),
                    (0.942378, (0.002463, 0.049449, 0.220096, 1.0)),
                    (1.0, (0.002405, 0.02925, 0.146109, 1.0)),
                ),
            )
        with g.Frame("Greens"):
            _string_8 = g.String(
                string="Greens: Sequential. 10 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_8 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.966503, 0.990422, 0.95535, 1.0)),
                    (0.173472, (0.705519, 0.880287, 0.652414, 1.0)),
                    (0.326861, (0.424202, 0.746674, 0.389969, 1.0)),
                    (0.46573, (0.209245, 0.587304, 0.19695, 1.0)),
                    (0.59334, (0.056124, 0.439663, 0.119367, 1.0)),
                    (0.68753, (0.031519, 0.32651, 0.081638, 1.0)),
                    (0.767024, (0.012005, 0.237971, 0.053013, 1.0)),
                    (0.850312, (0.0, 0.166249, 0.027706, 1.0)),
                    (0.902383, (0.0, 0.136296, 0.02254, 1.0)),
                    (1.0, (0.0, 0.041496, 0.00885, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("Oranges"):
            _string_9 = g.String(
                string="Oranges: Sequential. 11 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_9 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.998162, 0.943119, 0.873575, 1.0)),
                    (0.129285, (0.994972, 0.759168, 0.570782, 1.0)),
                    (0.204806, (0.980907, 0.718769, 0.474742, 1.0)),
                    (0.335093, (0.977088, 0.476131, 0.170929, 1.0)),
                    (0.472924, (1.0, 0.289435, 0.049796, 1.0)),
                    (0.590997, (0.919325, 0.161088, 0.006838, 1.0)),
                    (0.676042, (0.80304, 0.103794, 0.003496, 1.0)),
                    (0.747512, (0.702998, 0.061345, 0.0, 1.0)),
                    (0.800691, (0.547798, 0.05324, 0.000359, 1.0)),
                    (0.868865, (0.402982, 0.038399, 0.001144, 1.0)),
                    (1.0, (0.176855, 0.017092, 0.001087, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("Reds"):
            _string_10 = g.String(
                string="Reds: Sequential. 13 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_10 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.998959, 0.911461, 0.868681, 1.0)),
                    (0.121783, (0.992548, 0.748993, 0.647298, 1.0)),
                    (0.225271, (0.975111, 0.539402, 0.401558, 1.0)),
                    (0.311482, (0.97494, 0.38511, 0.250778, 1.0)),
                    (0.387279, (0.970297, 0.267799, 0.153396, 1.0)),
                    (0.486811, (0.969494, 0.15481, 0.074861, 1.0)),
                    (0.559143, (0.918792, 0.08691, 0.044058, 1.0)),
                    (0.62215, (0.86589, 0.044199, 0.025727, 1.0)),
                    (0.691798, (0.716226, 0.020882, 0.017407, 1.0)),
                    (0.750492, (0.591738, 0.008821, 0.012163, 1.0)),
                    (0.87451, (0.3736, 0.004809, 0.007438, 1.0)),
                    (0.942658, (0.226102, 0.002012, 0.005475, 1.0)),
                    (1.0, (0.134217, 0.000052, 0.003999, 1.0)),
                ),
            )
        with g.Frame("YlOrBr"):
            _string_11 = g.String(
                string="YlOrBr: Sequential. 13 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_11 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.999187, 1.0, 0.777707, 1.0)),
                    (0.132726, (0.999995, 0.92219, 0.482126, 1.0)),
                    (0.25098, (0.990341, 0.76396, 0.277625, 1.0)),
                    (0.3125, (0.992286, 0.655285, 0.159578, 1.0)),
                    (0.37121, (0.989224, 0.55707, 0.080231, 1.0)),
                    (0.439566, (0.993986, 0.422641, 0.043374, 1.0)),
                    (0.501023, (0.986543, 0.313863, 0.021427, 1.0)),
                    (0.601863, (0.870695, 0.184123, 0.008736, 1.0)),
                    (0.649634, (0.791772, 0.13902, 0.00523, 1.0)),
                    (0.746691, (0.606455, 0.072667, 0.000654, 1.0)),
                    (0.846165, (0.371521, 0.040553, 0.001034, 1.0)),
                    (0.92001, (0.237824, 0.027205, 0.001459, 1.0)),
                    (1.0, (0.13089, 0.018824, 0.001804, 1.0)),
                ),
            )
        with g.Frame("YlOrRd"):
            _string_12 = g.String(
                string="YlOrRd: Sequential. 14 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_12 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.999992, 1.0, 0.597977, 1.0)),
                    (0.123328, (1.0, 0.844676, 0.349616, 1.0)),
                    (0.237138, (0.991688, 0.713497, 0.192726, 1.0)),
                    (0.323335, (0.990584, 0.543101, 0.108901, 1.0)),
                    (0.378469, (0.991873, 0.435601, 0.070421, 1.0)),
                    (0.498114, (0.980933, 0.265895, 0.045518, 1.0)),
                    (0.570252, (0.97965, 0.140346, 0.031527, 1.0)),
                    (0.629039, (0.969335, 0.070919, 0.022518, 1.0)),
                    (0.695189, (0.852348, 0.029396, 0.016182, 1.0)),
                    (0.746788, (0.774378, 0.010596, 0.011643, 1.0)),
                    (0.800386, (0.654811, 0.004839, 0.014301, 1.0)),
                    (0.873309, (0.510363, 0.000043, 0.019332, 1.0)),
                    (0.939568, (0.336766, 0.000002, 0.0194, 1.0)),
                    (1.0, (0.214137, 0.0, 0.019421, 1.0)),
                ),
            )
        with g.Frame("OrRd"):
            _string_13 = g.String(
                string="OrRd: Sequential. 14 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_13 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.999998, 0.932983, 0.834433, 1.0)),
                    (0.158205, (0.98821, 0.768783, 0.505365, 1.0)),
                    (0.249215, (0.983211, 0.656918, 0.340619, 1.0)),
                    (0.374291, (0.980819, 0.495542, 0.22988, 1.0)),
                    (0.446472, (0.979572, 0.354694, 0.146535, 1.0)),
                    (0.504646, (0.969183, 0.256505, 0.096887, 1.0)),
                    (0.612751, (0.876607, 0.139317, 0.068985, 1.0)),
                    (0.691781, (0.764028, 0.064818, 0.031648, 1.0)),
                    (0.742942, (0.690084, 0.032333, 0.014975, 1.0)),
                    (0.786006, (0.609544, 0.015675, 0.007861, 1.0)),
                    (0.825757, (0.533123, 0.006157, 0.003666, 1.0)),
                    (0.873295, (0.452308, 0.000028, 0.000068, 1.0)),
                    (0.939556, (0.312695, 0.0, 0.0, 1.0)),
                    (1.0, (0.210967, 0.0, 0.0, 1.0)),
                ),
            )
        with g.Frame("PuRd"):
            _string_14 = g.String(
                string="PuRd: Sequential. 13 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_14 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.927956, 0.904336, 0.947535, 1.0)),
                    (0.121712, (0.80426, 0.754873, 0.864524, 1.0)),
                    (0.248532, (0.656216, 0.483818, 0.704446, 1.0)),
                    (0.380372, (0.585379, 0.283741, 0.561592, 1.0)),
                    (0.479168, (0.716592, 0.149937, 0.459269, 1.0)),
                    (0.563466, (0.774081, 0.060384, 0.338174, 1.0)),
                    (0.62154, (0.798792, 0.022616, 0.256739, 1.0)),
                    (0.688882, (0.705407, 0.012245, 0.158773, 1.0)),
                    (0.748, (0.619379, 0.006051, 0.09355, 1.0)),
                    (0.81375, (0.445764, 0.002594, 0.073129, 1.0)),
                    (0.874629, (0.312785, 0.000056, 0.055478, 1.0)),
                    (0.942647, (0.205123, 0.0, 0.028383, 1.0)),
                    (1.0, (0.13477, 0.0, 0.013395, 1.0)),
                ),
            )
        with g.Frame("RdPu"):
            _string_15 = g.String(
                string="RdPu: Sequential. 14 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_15 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.997729, 0.931405, 0.902237, 1.0)),
                    (0.250752, (0.971298, 0.5528, 0.529357, 1.0)),
                    (0.376471, (0.956959, 0.340569, 0.459171, 1.0)),
                    (0.444078, (0.942038, 0.216239, 0.402209, 1.0)),
                    (0.501663, (0.92819, 0.134761, 0.355053, 1.0)),
                    (0.566621, (0.814723, 0.071356, 0.33127, 1.0)),
                    (0.624677, (0.723459, 0.033614, 0.309363, 1.0)),
                    (0.671637, (0.59886, 0.014722, 0.268925, 1.0)),
                    (0.707933, (0.514949, 0.005845, 0.240068, 1.0)),
                    (0.749815, (0.42116, 0.000207, 0.20837, 1.0)),
                    (0.821201, (0.278281, 0.000366, 0.195119, 1.0)),
                    (0.883696, (0.181369, 0.000249, 0.182029, 1.0)),
                    (0.946807, (0.11033, 0.000155, 0.160153, 1.0)),
                    (1.0, (0.065962, 0.000001, 0.144354, 1.0)),
                ),
            )
        with g.Frame("BuPu"):
            _string_16 = g.String(
                string="BuPu: Sequential. 12 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.7/255 in sRGB."
            )
            color_ramp_16 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.928649, 0.972453, 0.981656, 1.0)),
                    (0.121001, (0.749204, 0.842762, 0.90756, 1.0)),
                    (0.249705, (0.518326, 0.650108, 0.79096, 1.0)),
                    (0.370186, (0.345451, 0.506304, 0.704403, 1.0)),
                    (0.495432, (0.263396, 0.307999, 0.569814, 1.0)),
                    (0.597283, (0.263126, 0.173885, 0.463353, 1.0)),
                    (0.680252, (0.256236, 0.096937, 0.394851, 1.0)),
                    (0.767411, (0.242718, 0.04123, 0.318542, 1.0)),
                    (0.825039, (0.229307, 0.016038, 0.248842, 1.0)),
                    (0.872039, (0.221137, 0.004881, 0.203281, 1.0)),
                    (0.942612, (0.128159, 0.002025, 0.119345, 1.0)),
                    (1.0, (0.073321, 0.000044, 0.069585, 1.0)),
                ),
            )
        with g.Frame("GnBu"):
            _string_17 = g.String(
                string="GnBu: Sequential. 12 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_17 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.928487, 0.974233, 0.87008, 1.0)),
                    (0.122831, (0.746898, 0.896191, 0.709667, 1.0)),
                    (0.24627, (0.606102, 0.834376, 0.561125, 1.0)),
                    (0.37459, (0.388792, 0.720708, 0.4614, 1.0)),
                    (0.467999, (0.237901, 0.636522, 0.527278, 1.0)),
                    (0.548178, (0.141211, 0.544545, 0.5901, 1.0)),
                    (0.626456, (0.073462, 0.445812, 0.650666, 1.0)),
                    (0.726375, (0.03029, 0.290863, 0.537597, 1.0)),
                    (0.806294, (0.010466, 0.199173, 0.465875, 1.0)),
                    (0.871904, (0.002383, 0.139862, 0.415964, 1.0)),
                    (0.942248, (0.002475, 0.084453, 0.298689, 1.0)),
                    (1.0, (0.002397, 0.05077, 0.218743, 1.0)),
                ),
            )
        with g.Frame("PuBu"):
            _string_18 = g.String(
                string="PuBu: Sequential. 14 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.6/255 in sRGB."
            )
            color_ramp_18 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.998893, 0.929372, 0.964389, 1.0)),
                    (0.12168, (0.841868, 0.801643, 0.889821, 1.0)),
                    (0.250782, (0.627062, 0.635679, 0.790719, 1.0)),
                    (0.325088, (0.470627, 0.556307, 0.739948, 1.0)),
                    (0.441521, (0.257873, 0.448294, 0.664227, 1.0)),
                    (0.54333, (0.111838, 0.354844, 0.59091, 1.0)),
                    (0.589026, (0.063643, 0.309129, 0.552922, 1.0)),
                    (0.62579, (0.035772, 0.278497, 0.527243, 1.0)),
                    (0.671585, (0.016982, 0.230145, 0.490121, 1.0)),
                    (0.708836, (0.007247, 0.197376, 0.465222, 1.0)),
                    (0.747875, (0.001541, 0.162115, 0.433134, 1.0)),
                    (0.875578, (0.001249, 0.100926, 0.262838, 1.0)),
                    (0.942802, (0.00086, 0.063408, 0.16107, 1.0)),
                    (1.0, (0.000623, 0.039195, 0.096612, 1.0)),
                ),
            )
        with g.Frame("YlGnBu"):
            _string_19 = g.String(
                string="YlGnBu: Sequential. 15 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_19 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (1.0, 0.999615, 0.689003, 1.0)),
                    (0.122407, (0.847975, 0.940256, 0.440319, 1.0)),
                    (0.248663, (0.569798, 0.816092, 0.455483, 1.0)),
                    (0.32286, (0.333016, 0.690504, 0.479973, 1.0)),
                    (0.381807, (0.196977, 0.600972, 0.499486, 1.0)),
                    (0.447761, (0.102507, 0.525094, 0.528267, 1.0)),
                    (0.499412, (0.052247, 0.467468, 0.55225, 1.0)),
                    (0.566038, (0.026747, 0.362731, 0.538083, 1.0)),
                    (0.622462, (0.012388, 0.284952, 0.528705, 1.0)),
                    (0.69161, (0.014066, 0.17892, 0.452202, 1.0)),
                    (0.750683, (0.016075, 0.110054, 0.390302, 1.0)),
                    (0.816832, (0.017367, 0.063116, 0.33901, 1.0)),
                    (0.873072, (0.018324, 0.034515, 0.296507, 1.0)),
                    (0.942723, (0.007337, 0.020521, 0.170706, 1.0)),
                    (1.0, (0.002284, 0.012167, 0.096362, 1.0)),
                ),
            )
        with g.Frame("PuBuGn"):
            _string_20 = g.String(
                string="PuBuGn: Sequential. 16 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_20 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.99878, 0.928779, 0.964299, 1.0)),
                    (0.122791, (0.8404, 0.762045, 0.87254, 1.0)),
                    (0.242593, (0.639845, 0.64392, 0.795766, 1.0)),
                    (0.376471, (0.373923, 0.506413, 0.707072, 1.0)),
                    (0.444058, (0.225182, 0.445218, 0.660863, 1.0)),
                    (0.498952, (0.135667, 0.397249, 0.624629, 1.0)),
                    (0.549672, (0.086914, 0.347033, 0.584005, 1.0)),
                    (0.588935, (0.057892, 0.31025, 0.554561, 1.0)),
                    (0.624929, (0.036267, 0.278944, 0.52634, 1.0)),
                    (0.668433, (0.017042, 0.257255, 0.417998, 1.0)),
                    (0.708064, (0.006379, 0.238628, 0.332877, 1.0)),
                    (0.749151, (0.000569, 0.219624, 0.254007, 1.0)),
                    (0.812644, (0.000494, 0.182679, 0.16507, 1.0)),
                    (0.873115, (0.000281, 0.150395, 0.100345, 1.0)),
                    (0.941804, (0.000327, 0.09607, 0.061079, 1.0)),
                    (1.0, (0.000304, 0.060749, 0.036519, 1.0)),
                ),
            )
        with g.Frame("BuGn"):
            _string_21 = g.String(
                string="BuGn: Sequential. 11 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_21 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.954641, 0.981569, 0.980144, 1.0)),
                    (0.126833, (0.768722, 0.911847, 0.965415, 1.0)),
                    (0.241495, (0.626441, 0.846112, 0.791853, 1.0)),
                    (0.32362, (0.423658, 0.758462, 0.67377, 1.0)),
                    (0.43637, (0.192556, 0.598984, 0.473539, 1.0)),
                    (0.542937, (0.0909, 0.499969, 0.291875, 1.0)),
                    (0.641224, (0.044147, 0.408414, 0.156082, 1.0)),
                    (0.734948, (0.018039, 0.264838, 0.057463, 1.0)),
                    (0.841452, (0.0, 0.171639, 0.031049, 1.0)),
                    (0.904419, (0.0, 0.133266, 0.020828, 1.0)),
                    (1.0, (0.0, 0.042298, 0.009413, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("YlGn"):
            _string_22 = g.String(
                string="YlGn: Sequential. 13 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_22 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.997362, 0.9993, 0.777902, 1.0)),
                    (0.119603, (0.936557, 0.975556, 0.491467, 1.0)),
                    (0.242017, (0.703985, 0.88006, 0.37072, 1.0)),
                    (0.412798, (0.332687, 0.674789, 0.24293, 1.0)),
                    (0.499337, (0.185945, 0.563556, 0.191883, 1.0)),
                    (0.569366, (0.099328, 0.474209, 0.141778, 1.0)),
                    (0.6244, (0.052367, 0.407293, 0.109781, 1.0)),
                    (0.69368, (0.029828, 0.301224, 0.076751, 1.0)),
                    (0.757794, (0.014651, 0.222372, 0.054062, 1.0)),
                    (0.812427, (0.005579, 0.181607, 0.047101, 1.0)),
                    (0.873533, (0.00001, 0.13875, 0.038103, 1.0)),
                    (0.92889, (0.000026, 0.098758, 0.030667, 1.0)),
                    (1.0, (0.0, 0.058852, 0.022041, 1.0)),
                ),
            )
        with g.Frame("binary"):
            _string_23 = g.String(
                string="binary: Sequential (2). 9 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_23 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.993647, 0.993647, 0.993647, 1.0)),
                    (0.171583, (0.647458, 0.647458, 0.647458, 1.0)),
                    (0.337154, (0.392005, 0.392005, 0.392005, 1.0)),
                    (0.49659, (0.213342, 0.213342, 0.213342, 1.0)),
                    (0.639992, (0.104056, 0.104056, 0.104056, 1.0)),
                    (0.757903, (0.046354, 0.046354, 0.046354, 1.0)),
                    (0.851687, (0.018494, 0.018494, 0.018494, 1.0)),
                    (0.924093, (0.006271, 0.006271, 0.006271, 1.0)),
                    (1.0, (0.0, 0.0, 0.0, 1.0)),
                ),
            )
        with g.Frame("gray"):
            _string_24 = g.String(
                string="gray: Sequential (2). 9 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_24 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.0, 1.0)),
                    (0.075907, (0.006271, 0.006271, 0.006271, 1.0)),
                    (0.148313, (0.018494, 0.018494, 0.018494, 1.0)),
                    (0.242097, (0.046354, 0.046354, 0.046354, 1.0)),
                    (0.360008, (0.104056, 0.104056, 0.104056, 1.0)),
                    (0.50341, (0.213342, 0.213342, 0.213342, 1.0)),
                    (0.662846, (0.392005, 0.392005, 0.392005, 1.0)),
                    (0.828417, (0.647458, 0.647458, 0.647458, 1.0)),
                    (1.0, (0.993647, 0.993647, 0.993647, 1.0)),
                ),
            )
        with g.Frame("bone"):
            _string_25 = g.String(
                string="bone: Sequential (2). 13 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_25 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.000044, 0.000044, 0.0, 1.0)),
                    (0.05316, (0.003538, 0.003538, 0.005176, 1.0)),
                    (0.097272, (0.007676, 0.007676, 0.012729, 1.0)),
                    (0.155336, (0.016182, 0.016182, 0.029071, 1.0)),
                    (0.224373, (0.031501, 0.031501, 0.059623, 1.0)),
                    (0.29521, (0.05375, 0.053753, 0.10505, 1.0)),
                    (0.367171, (0.083641, 0.083596, 0.16688, 1.0)),
                    (0.442309, (0.123152, 0.140198, 0.224473, 1.0)),
                    (0.532334, (0.182633, 0.23136, 0.306622, 1.0)),
                    (0.636765, (0.268923, 0.372215, 0.421288, 1.0)),
                    (0.751326, (0.388943, 0.572735, 0.572424, 1.0)),
                    (0.869259, (0.634681, 0.756721, 0.756763, 1.0)),
                    (1.0, (0.993028, 0.997225, 0.997208, 1.0)),
                ),
            )
        with g.Frame("pink"):
            _string_26 = g.String(
                string="pink: Sequential (2). 8 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_26 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.012843, 0.00004, 0.000039, 1.0)),
                    (0.009317, (0.022691, 0.004683, 0.004686, 1.0)),
                    (0.015432, (0.030323, 0.010016, 0.010002, 1.0)),
                    (0.131004, (0.180013, 0.07016, 0.070235, 1.0)),
                    (0.366382, (0.535124, 0.207731, 0.207389, 1.0)),
                    (0.534233, (0.658477, 0.460624, 0.313624, 1.0)),
                    (0.749797, (0.811297, 0.814248, 0.459569, 1.0)),
                    (1.0, (1.0, 0.998611, 0.995368, 1.0)),
                ),
            )
        with g.Frame("spring"):
            _string_27 = g.String(
                string="spring: Sequential (2). 10 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_27 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (1.0, 0.0, 1.0, 1.0)),
                    (0.028723, (1.0, 0.0014, 0.953548, 1.0)),
                    (0.118177, (1.0, 0.009557, 0.747718, 1.0)),
                    (0.244737, (1.0, 0.041422, 0.529613, 1.0)),
                    (0.407004, (1.0, 0.12398, 0.298246, 1.0)),
                    (0.592995, (1.0, 0.298246, 0.12398, 1.0)),
                    (0.755263, (1.0, 0.529613, 0.041422, 1.0)),
                    (0.881823, (1.0, 0.747718, 0.009557, 1.0)),
                    (0.971277, (1.0, 0.953548, 0.0014, 1.0)),
                    (1.0, (1.0, 1.0, 0.0, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("summer"):
            _string_28 = g.String(
                string="summer: Sequential (2). 9 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_28 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.213833, 0.132916, 1.0)),
                    (0.075907, (0.006271, 0.250619, 0.132873, 1.0)),
                    (0.148313, (0.018494, 0.288904, 0.13285, 1.0)),
                    (0.242097, (0.046354, 0.343164, 0.132874, 1.0)),
                    (0.360008, (0.104056, 0.419214, 0.132871, 1.0)),
                    (0.502922, (0.21297, 0.523719, 0.132867, 1.0)),
                    (0.660029, (0.388319, 0.654637, 0.13287, 1.0)),
                    (0.82781, (0.646232, 0.813715, 0.132867, 1.0)),
                    (1.0, (0.993647, 0.998447, 0.132869, 1.0)),
                ),
            )
        with g.Frame("autumn"):
            _string_29 = g.String(
                string="autumn: Sequential (2). 9 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_29 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.999999, 0.0, 0.0, 1.0)),
                    (0.07635, (1.0, 0.006308, 0.0, 1.0)),
                    (0.152064, (1.0, 0.019277, 0.0, 1.0)),
                    (0.249817, (1.0, 0.049263, 0.0, 1.0)),
                    (0.370739, (1.0, 0.110729, 0.0, 1.0)),
                    (0.511431, (1.0, 0.221085, 0.0, 1.0)),
                    (0.6674, (1.0, 0.398198, 0.0, 1.0)),
                    (0.831841, (1.0, 0.653657, 0.0, 1.0)),
                    (1.0, (1.0, 0.993955, 0.0, 1.0)),
                ),
            )
        with g.Frame("winter"):
            _string_30 = g.String(
                string="winter: Sequential (2). 9 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_30 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.999715, 1.0)),
                    (0.075907, (0.0, 0.006271, 0.915532, 1.0)),
                    (0.148313, (0.0, 0.018494, 0.839154, 1.0)),
                    (0.242097, (0.0, 0.046354, 0.745816, 1.0)),
                    (0.360008, (0.0, 0.104056, 0.637381, 1.0)),
                    (0.50341, (0.0, 0.213342, 0.51869, 1.0)),
                    (0.662846, (0.0, 0.392005, 0.40335, 1.0)),
                    (0.828418, (0.0, 0.647459, 0.300886, 1.0)),
                    (1.0, (0.0, 0.993647, 0.212958, 1.0)),
                ),
            )
        with g.Frame("cool"):
            _string_31 = g.String(
                string="cool: Sequential (2). 11 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_31 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 1.0, 1.0, 1.0)),
                    (0.02798, (0.00144, 0.954881, 1.0, 1.0)),
                    (0.106888, (0.008665, 0.767282, 1.0, 1.0)),
                    (0.210214, (0.031755, 0.588316, 1.0, 1.0)),
                    (0.343094, (0.087228, 0.380671, 1.0, 1.0)),
                    (0.498481, (0.204164, 0.207138, 1.0, 1.0)),
                    (0.654671, (0.37741, 0.088382, 1.0, 1.0)),
                    (0.78899, (0.587036, 0.031875, 1.0, 1.0)),
                    (0.893127, (0.767585, 0.008663, 1.0, 1.0)),
                    (0.972021, (0.954766, 0.00144, 1.0, 1.0)),
                    (1.0, (1.0, 0.0, 1.0, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("Wistia"):
            _string_32 = g.String(
                string="Wistia: Sequential (2). 9 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.7/255 in sRGB."
            )
            color_ramp_32 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.775629, 0.999837, 0.192224, 1.0)),
                    (0.080033, (0.843738, 0.935376, 0.103309, 1.0)),
                    (0.149904, (0.905974, 0.88114, 0.050811, 1.0)),
                    (0.205104, (0.956937, 0.839506, 0.023773, 1.0)),
                    (0.249076, (0.999574, 0.807613, 0.01016, 1.0)),
                    (0.357841, (1.0, 0.6669, 0.004555, 1.0)),
                    (0.501812, (1.0, 0.504841, 0.0, 1.0)),
                    (0.809211, (0.996583, 0.313008, 0.0, 1.0)),
                    (1.0, (0.971423, 0.210817, 0.0, 1.0)),
                ),
            )
        with g.Frame("hot"):
            _string_33 = g.String(
                string="hot: Sequential (2). 24 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_33 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.0, 1.0)),
                    (0.052744, (0.020264, 0.0, 0.0, 1.0)),
                    (0.11564, (0.085942, 0.0, 0.0, 1.0)),
                    (0.179229, (0.220483, 0.0, 0.0, 1.0)),
                    (0.234334, (0.374962, 0.0, 0.0, 1.0)),
                    (0.285871, (0.598127, 0.0, 0.0, 1.0)),
                    (0.324016, (0.752955, 0.0, 0.0, 1.0)),
                    (0.359515, (0.996626, 0.0, 0.0, 1.0)),
                    (0.370502, (1.0, 0.0, 0.0, 1.0)),
                    (0.409718, (1.0, 0.008969, 0.0, 1.0)),
                    (0.463985, (1.0, 0.045103, 0.0, 1.0)),
                    (0.528851, (1.0, 0.145771, 0.0, 1.0)),
                    (0.595655, (1.0, 0.307843, 0.0, 1.0)),
                    (0.656842, (1.0, 0.554097, 0.0, 1.0)),
                    (0.700709, (1.0, 0.730032, 0.0, 1.0)),
                    (0.739728, (1.0, 0.99518, 0.0, 1.0)),
                    (0.753282, (1.0, 1.0, 0.0, 1.0)),
                    (0.786058, (1.0, 1.0, 0.014864, 1.0)),
                    (0.828688, (1.0, 1.0, 0.073427, 1.0)),
                    (0.878484, (1.0, 1.0, 0.217872, 1.0)),
                    (0.924135, (1.0, 1.0, 0.449038, 1.0)),
                    (0.960037, (1.0, 1.0, 0.660379, 1.0)),
                    (0.99134, (1.0, 1.0, 0.950063, 1.0)),
                    (1.0, (1.0, 1.0, 1.0, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("afmhot"):
            _string_34 = g.String(
                string="afmhot: Sequential (2). 19 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_34 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.0, 1.0)),
                    (0.034757, (0.005574, 0.0, 0.0, 1.0)),
                    (0.073142, (0.017862, 0.0, 0.0, 1.0)),
                    (0.123507, (0.048092, 0.0, 0.0, 1.0)),
                    (0.182464, (0.107213, 0.0, 0.0, 1.0)),
                    (0.249013, (0.209079, 0.0, 0.0, 1.0)),
                    (0.285045, (0.286068, 0.005571, 0.0, 1.0)),
                    (0.324569, (0.376622, 0.018447, 0.0, 1.0)),
                    (0.377785, (0.528847, 0.051375, 0.0, 1.0)),
                    (0.437023, (0.733896, 0.113084, 0.0, 1.0)),
                    (0.500767, (0.999581, 0.212661, 0.0, 1.0)),
                    (0.53915, (1.0, 0.294917, 0.006542, 1.0)),
                    (0.575327, (1.0, 0.37895, 0.018999, 1.0)),
                    (0.622181, (1.0, 0.511497, 0.047195, 1.0)),
                    (0.680851, (1.0, 0.710562, 0.105122, 1.0)),
                    (0.751287, (1.0, 1.0, 0.212555, 1.0)),
                    (0.832649, (1.0, 1.0, 0.394858, 1.0)),
                    (0.918058, (1.0, 1.0, 0.661065, 1.0)),
                    (1.0, (1.0, 1.0, 0.994688, 1.0)),
                ),
            )
        with g.Frame("gist_heat"):
            _string_35 = g.String(
                string="gist_heat: Sequential (2). 24 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_35 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.0, 1.0)),
                    (0.041854, (0.004985, 0.0, 0.0, 1.0)),
                    (0.07737, (0.012301, 0.0, 0.0, 1.0)),
                    (0.124213, (0.028252, 0.0, 0.0, 1.0)),
                    (0.182887, (0.059901, 0.0, 0.0, 1.0)),
                    (0.253371, (0.117346, 0.0, 0.0, 1.0)),
                    (0.334159, (0.212616, 0.0, 0.0, 1.0)),
                    (0.416815, (0.345985, 0.0, 0.0, 1.0)),
                    (0.499599, (0.51865, 0.0, 0.0, 1.0)),
                    (0.532205, (0.602133, 0.005091, 0.0, 1.0)),
                    (0.564103, (0.685615, 0.014334, 0.0, 1.0)),
                    (0.609898, (0.813482, 0.038278, 0.0, 1.0)),
                    (0.664155, (0.995234, 0.086214, 0.0, 1.0)),
                    (0.709525, (1.0, 0.146002, 0.0, 1.0)),
                    (0.750761, (1.0, 0.214079, 0.0, 1.0)),
                    (0.770301, (1.0, 0.254298, 0.006862, 1.0)),
                    (0.787665, (1.0, 0.289818, 0.018913, 1.0)),
                    (0.813316, (1.0, 0.350104, 0.050425, 1.0)),
                    (0.841387, (1.0, 0.423104, 0.107777, 1.0)),
                    (0.8717, (1.0, 0.511609, 0.19951, 1.0)),
                    (0.899811, (1.0, 0.602528, 0.315315, 1.0)),
                    (0.928336, (1.0, 0.70353, 0.4641, 1.0)),
                    (0.963003, (1.0, 0.838631, 0.691072, 1.0)),
                    (1.0, (1.0, 0.999015, 0.995937, 1.0)),
                ),
            )
        with g.Frame("copper"):
            _string_36 = g.String(
                string="copper: Sequential (2). 11 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_36 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.000016, 0.000069, 0.000195, 1.0)),
                    (0.054119, (0.005372, 0.003201, 0.002053, 1.0)),
                    (0.10514, (0.014756, 0.007278, 0.004077, 1.0)),
                    (0.171655, (0.035971, 0.015811, 0.007776, 1.0)),
                    (0.252475, (0.077686, 0.031715, 0.014277, 1.0)),
                    (0.340706, (0.14601, 0.056907, 0.024145, 1.0)),
                    (0.439774, (0.253342, 0.095673, 0.03891, 1.0)),
                    (0.555569, (0.424598, 0.15657, 0.061607, 1.0)),
                    (0.677331, (0.663301, 0.240548, 0.09237, 1.0)),
                    (0.812466, (1.0, 0.357282, 0.134703, 1.0)),
                    (1.0, (1.0, 0.567996, 0.210206, 1.0)),
                ),
            )
        with g.Frame("PiYG"):
            _string_37 = g.String(
                string="PiYG: Diverging. 18 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_37 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.267574, 0.000171, 0.08296, 1.0)),
                    (0.067056, (0.446083, 0.005963, 0.156295, 1.0)),
                    (0.103687, (0.567642, 0.011956, 0.211809, 1.0)),
                    (0.134226, (0.611533, 0.041375, 0.268705, 1.0)),
                    (0.16353, (0.668172, 0.089886, 0.331333, 1.0)),
                    (0.252985, (0.811063, 0.314245, 0.5577, 1.0)),
                    (0.298039, (0.875099, 0.457572, 0.695844, 1.0)),
                    (0.402947, (0.98308, 0.748041, 0.866549, 1.0)),
                    (0.49915, (0.929457, 0.92851, 0.927558, 1.0)),
                    (0.599794, (0.789468, 0.912358, 0.626583, 1.0)),
                    (0.655969, (0.605058, 0.823058, 0.379909, 1.0)),
                    (0.704098, (0.464338, 0.743223, 0.22408, 1.0)),
                    (0.754536, (0.315815, 0.607167, 0.115436, 1.0)),
                    (0.798364, (0.213347, 0.505366, 0.052933, 1.0)),
                    (0.852343, (0.128095, 0.38022, 0.029354, 1.0)),
                    (0.899802, (0.073307, 0.286062, 0.014955, 1.0)),
                    (0.954513, (0.039198, 0.189492, 0.01217, 1.0)),
                    (1.0, (0.019936, 0.126719, 0.009631, 1.0)),
                ),
            )
        with g.Frame("PRGn"):
            _string_38 = g.String(
                string="PRGn: Diverging. 22 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.7/255 in sRGB."
            )
            color_ramp_38 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.050814, 0.000056, 0.069838, 1.0)),
                    (0.032991, (0.083165, 0.004182, 0.109712, 1.0)),
                    (0.062084, (0.120087, 0.01008, 0.154289, 1.0)),
                    (0.100191, (0.180327, 0.02266, 0.226114, 1.0)),
                    (0.128171, (0.215754, 0.046748, 0.270614, 1.0)),
                    (0.1628, (0.26139, 0.091391, 0.333883, 1.0)),
                    (0.202372, (0.321757, 0.164935, 0.409607, 1.0)),
                    (0.248821, (0.41771, 0.251875, 0.509014, 1.0)),
                    (0.319833, (0.584125, 0.421628, 0.661774, 1.0)),
                    (0.402134, (0.801928, 0.660088, 0.806856, 1.0)),
                    (0.499794, (0.927477, 0.926398, 0.926317, 1.0)),
                    (0.6, (0.691448, 0.873632, 0.648013, 1.0)),
                    (0.652857, (0.513192, 0.780997, 0.477822, 1.0)),
                    (0.699202, (0.381892, 0.710526, 0.352146, 1.0)),
                    (0.734902, (0.258304, 0.598221, 0.253004, 1.0)),
                    (0.767623, (0.168509, 0.506025, 0.177954, 1.0)),
                    (0.800834, (0.099508, 0.420211, 0.117495, 1.0)),
                    (0.836454, (0.055345, 0.32327, 0.083629, 1.0)),
                    (0.866425, (0.028933, 0.255211, 0.05936, 1.0)),
                    (0.89924, (0.010679, 0.187338, 0.038136, 1.0)),
                    (0.951283, (0.003899, 0.108543, 0.021437, 1.0)),
                    (1.0, (0.00008, 0.056857, 0.010755, 1.0)),
                ),
            )
        with g.Frame("BrBG"):
            _string_39 = g.String(
                string="BrBG: Diverging. 22 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_39 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.087617, 0.029262, 0.001456, 1.0)),
                    (0.046023, (0.153236, 0.049396, 0.002311, 1.0)),
                    (0.100716, (0.261248, 0.082119, 0.002901, 1.0)),
                    (0.143828, (0.360985, 0.131521, 0.009476, 1.0)),
                    (0.199573, (0.517652, 0.216744, 0.025168, 1.0)),
                    (0.229557, (0.581974, 0.296246, 0.057338, 1.0)),
                    (0.264362, (0.654035, 0.404594, 0.116233, 1.0)),
                    (0.300464, (0.740581, 0.539588, 0.204092, 1.0)),
                    (0.343752, (0.812029, 0.647478, 0.327322, 1.0)),
                    (0.400072, (0.926356, 0.805799, 0.539928, 1.0)),
                    (0.500399, (0.905835, 0.914444, 0.905652, 1.0)),
                    (0.6, (0.565415, 0.819967, 0.7868, 1.0)),
                    (0.65396, (0.351479, 0.706952, 0.638781, 1.0)),
                    (0.705259, (0.198194, 0.593911, 0.517802, 1.0)),
                    (0.75467, (0.092461, 0.427691, 0.37712, 1.0)),
                    (0.797594, (0.036295, 0.315577, 0.27952, 1.0)),
                    (0.832829, (0.017122, 0.240587, 0.21097, 1.0)),
                    (0.861538, (0.007203, 0.189904, 0.163683, 1.0)),
                    (0.88757, (0.002086, 0.150681, 0.128376, 1.0)),
                    (0.902632, (0.00016, 0.128666, 0.107601, 1.0)),
                    (0.954502, (0.000223, 0.077441, 0.058269, 1.0)),
                    (1.0, (0.0, 0.044671, 0.029, 1.0)),
                ),
            )
        with g.Frame("PuOr"):
            _string_40 = g.String(
                string="PuOr: Diverging. 24 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.7/255 in sRGB."
            )
            color_ramp_40 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.211181, 0.043473, 0.002423, 1.0)),
                    (0.04587, (0.30684, 0.064919, 0.002158, 1.0)),
                    (0.10012, (0.448745, 0.096935, 0.001807, 1.0)),
                    (0.147839, (0.579747, 0.149073, 0.003894, 1.0)),
                    (0.201397, (0.748083, 0.224225, 0.006947, 1.0)),
                    (0.22519, (0.801699, 0.278052, 0.020439, 1.0)),
                    (0.248799, (0.855667, 0.332921, 0.042327, 1.0)),
                    (0.272899, (0.914874, 0.398705, 0.07442, 1.0)),
                    (0.300889, (0.983276, 0.480897, 0.125113, 1.0)),
                    (0.334946, (0.984694, 0.563922, 0.213538, 1.0)),
                    (0.369552, (0.989111, 0.656693, 0.333705, 1.0)),
                    (0.4, (0.989959, 0.744217, 0.465886, 1.0)),
                    (0.445824, (0.964595, 0.828563, 0.65344, 1.0)),
                    (0.500331, (0.927269, 0.927659, 0.927376, 1.0)),
                    (0.597328, (0.688818, 0.703961, 0.831302, 1.0)),
                    (0.683456, (0.481248, 0.446392, 0.677358, 1.0)),
                    (0.755778, (0.305877, 0.259833, 0.510526, 1.0)),
                    (0.794175, (0.224904, 0.181032, 0.422038, 1.0)),
                    (0.840658, (0.155098, 0.08686, 0.339294, 1.0)),
                    (0.875198, (0.114242, 0.041111, 0.28234, 1.0)),
                    (0.899638, (0.088514, 0.020019, 0.245691, 1.0)),
                    (0.934978, (0.061511, 0.009675, 0.169347, 1.0)),
                    (0.963866, (0.043441, 0.004306, 0.11841, 1.0)),
                    (1.0, (0.025976, 0.000038, 0.069597, 1.0)),
                ),
            )
        with g.Frame("RdGy"):
            _string_41 = g.String(
                string="RdGy: Diverging. 16 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_41 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.13372, 0.0001, 0.013692, 1.0)),
                    (0.045532, (0.246855, 0.003165, 0.017944, 1.0)),
                    (0.101183, (0.444025, 0.008983, 0.024288, 1.0)),
                    (0.134598, (0.519772, 0.029307, 0.037872, 1.0)),
                    (0.169823, (0.596213, 0.067595, 0.054948, 1.0)),
                    (0.205758, (0.685758, 0.125601, 0.078927, 1.0)),
                    (0.249495, (0.781102, 0.220939, 0.134839, 1.0)),
                    (0.299812, (0.902426, 0.372045, 0.220754, 1.0)),
                    (0.342694, (0.94026, 0.505465, 0.345995, 1.0)),
                    (0.414885, (0.988474, 0.752658, 0.622388, 1.0)),
                    (0.501377, (0.996615, 0.996557, 0.996635, 1.0)),
                    (0.775189, (0.29216, 0.292175, 0.292155, 1.0)),
                    (0.846072, (0.1483, 0.1483, 0.1483, 1.0)),
                    (0.900479, (0.071582, 0.071582, 0.071582, 1.0)),
                    (0.957807, (0.027798, 0.027798, 0.027799, 1.0)),
                    (1.0, (0.009906, 0.009906, 0.009906, 1.0)),
                ),
            )
        with g.Frame("RdBu"):
            _string_42 = g.String(
                string="RdBu: Diverging. 21 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_42 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.133719, 0.000105, 0.013739, 1.0)),
                    (0.045555, (0.246915, 0.003158, 0.017902, 1.0)),
                    (0.100755, (0.442569, 0.008974, 0.024347, 1.0)),
                    (0.120245, (0.489625, 0.019764, 0.031384, 1.0)),
                    (0.143694, (0.535636, 0.037897, 0.04254, 1.0)),
                    (0.175248, (0.612355, 0.075331, 0.057691, 1.0)),
                    (0.208751, (0.692023, 0.131269, 0.081881, 1.0)),
                    (0.252756, (0.788659, 0.229466, 0.139981, 1.0)),
                    (0.300563, (0.906696, 0.374989, 0.22231, 1.0)),
                    (0.346946, (0.938198, 0.514814, 0.360042, 1.0)),
                    (0.399771, (0.985795, 0.705286, 0.56474, 1.0)),
                    (0.500932, (0.924359, 0.930147, 0.924957, 1.0)),
                    (0.6, (0.632958, 0.779645, 0.873395, 1.0)),
                    (0.652567, (0.433399, 0.664311, 0.797688, 1.0)),
                    (0.714372, (0.238375, 0.518175, 0.705179, 1.0)),
                    (0.758246, (0.124598, 0.387145, 0.616838, 1.0)),
                    (0.798331, (0.056381, 0.294364, 0.547004, 1.0)),
                    (0.848238, (0.03204, 0.205265, 0.480591, 1.0)),
                    (0.898986, (0.014992, 0.132017, 0.41061, 1.0)),
                    (0.954508, (0.005656, 0.064177, 0.223583, 1.0)),
                    (1.0, (0.00147, 0.028825, 0.117679, 1.0)),
                ),
            )
        with g.Frame("RdYlBu"):
            _string_43 = g.String(
                string="RdYlBu: Diverging. 22 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_43 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.375643, 0.0, 0.019414, 1.0)),
                    (0.033718, (0.466067, 0.005045, 0.019637, 1.0)),
                    (0.061881, (0.550651, 0.012362, 0.020015, 1.0)),
                    (0.100124, (0.67882, 0.028848, 0.020175, 1.0)),
                    (0.134856, (0.752625, 0.058937, 0.030239, 1.0)),
                    (0.169402, (0.833199, 0.102084, 0.042749, 1.0)),
                    (0.203022, (0.909201, 0.157053, 0.057383, 1.0)),
                    (0.2491, (0.940791, 0.263218, 0.083183, 1.0)),
                    (0.301061, (0.983451, 0.423251, 0.120185, 1.0)),
                    (0.343186, (0.984967, 0.551546, 0.177165, 1.0)),
                    (0.410293, (0.993645, 0.773445, 0.29581, 1.0)),
                    (0.499672, (0.99748, 0.994297, 0.514378, 1.0)),
                    (0.563574, (0.835729, 0.937168, 0.766863, 1.0)),
                    (0.598845, (0.742849, 0.892951, 0.936242, 1.0)),
                    (0.685279, (0.444652, 0.72728, 0.835905, 1.0)),
                    (0.752433, (0.268937, 0.542509, 0.72132, 1.0)),
                    (0.801046, (0.171314, 0.412015, 0.633855, 1.0)),
                    (0.852508, (0.104252, 0.275168, 0.538219, 1.0)),
                    (0.896805, (0.061193, 0.182693, 0.461068, 1.0)),
                    (0.924941, (0.050924, 0.129923, 0.413962, 1.0)),
                    (0.960514, (0.040792, 0.076508, 0.357147, 1.0)),
                    (1.0, (0.030573, 0.036035, 0.300214, 1.0)),
                ),
            )
        with g.Frame("RdYlGn"):
            _string_44 = g.String(
                string="RdYlGn: Diverging. 26 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.6/255 in sRGB."
            )
            color_ramp_44 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.375642, 0.0, 0.019413, 1.0)),
                    (0.033718, (0.466069, 0.005045, 0.019636, 1.0)),
                    (0.061881, (0.550594, 0.012362, 0.020018, 1.0)),
                    (0.099768, (0.677795, 0.028694, 0.020168, 1.0)),
                    (0.131547, (0.747152, 0.055638, 0.029172, 1.0)),
                    (0.165804, (0.82234, 0.096902, 0.041336, 1.0)),
                    (0.198303, (0.901851, 0.148823, 0.05514, 1.0)),
                    (0.226498, (0.925825, 0.20835, 0.070151, 1.0)),
                    (0.261082, (0.951052, 0.297572, 0.091155, 1.0)),
                    (0.301754, (0.983041, 0.426816, 0.120796, 1.0)),
                    (0.347747, (0.98623, 0.561919, 0.176918, 1.0)),
                    (0.400715, (0.990702, 0.745925, 0.2581, 1.0)),
                    (0.446364, (0.996035, 0.857238, 0.364972, 1.0)),
                    (0.5, (0.998741, 0.999124, 0.518391, 1.0)),
                    (0.553636, (0.826483, 0.92502, 0.36497, 1.0)),
                    (0.598925, (0.695831, 0.864938, 0.258945, 1.0)),
                    (0.648721, (0.526179, 0.777071, 0.197058, 1.0)),
                    (0.698446, (0.383207, 0.697076, 0.145063, 1.0)),
                    (0.747021, (0.243493, 0.601148, 0.134033, 1.0)),
                    (0.789176, (0.151476, 0.528954, 0.127944, 1.0)),
                    (0.825456, (0.085386, 0.455448, 0.11332, 1.0)),
                    (0.850043, (0.050113, 0.403416, 0.100408, 1.0)),
                    (0.879075, (0.02238, 0.351101, 0.088954, 1.0)),
                    (0.899892, (0.010018, 0.312677, 0.079932, 1.0)),
                    (0.950168, (0.003844, 0.213945, 0.056544, 1.0)),
                    (1.0, (0.000071, 0.137383, 0.037972, 1.0)),
                ),
            )
        with g.Frame("Spectral"):
            _string_45 = g.String(
                string="Spectral: Diverging. 25 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.7/255 in sRGB."
            )
            color_ramp_45 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.341614, 0.000319, 0.054527, 1.0)),
                    (0.022358, (0.403013, 0.004479, 0.059251, 1.0)),
                    (0.039812, (0.456032, 0.009718, 0.063457, 1.0)),
                    (0.057047, (0.510264, 0.0172, 0.067159, 1.0)),
                    (0.076343, (0.578176, 0.028837, 0.072449, 1.0)),
                    (0.101816, (0.670229, 0.048943, 0.078112, 1.0)),
                    (0.144619, (0.765684, 0.085195, 0.067519, 1.0)),
                    (0.201073, (0.906396, 0.152562, 0.056184, 1.0)),
                    (0.249127, (0.942, 0.263281, 0.083166, 1.0)),
                    (0.302782, (0.983491, 0.428607, 0.12152, 1.0)),
                    (0.351256, (0.986319, 0.57341, 0.181831, 1.0)),
                    (0.401271, (0.990575, 0.748163, 0.259161, 1.0)),
                    (0.449195, (0.996614, 0.863946, 0.372289, 1.0)),
                    (0.499499, (0.998195, 0.999091, 0.517312, 1.0)),
                    (0.598929, (0.793179, 0.91462, 0.312018, 1.0)),
                    (0.649408, (0.580683, 0.816142, 0.343342, 1.0)),
                    (0.711168, (0.365721, 0.701322, 0.374296, 1.0)),
                    (0.754901, (0.232431, 0.617813, 0.372446, 1.0)),
                    (0.799071, (0.132273, 0.539941, 0.376086, 1.0)),
                    (0.852047, (0.068862, 0.367938, 0.441183, 1.0)),
                    (0.895814, (0.033718, 0.254463, 0.503952, 1.0)),
                    (0.902311, (0.032662, 0.240512, 0.506271, 1.0)),
                    (0.935402, (0.053319, 0.172977, 0.452322, 1.0)),
                    (0.970271, (0.081567, 0.115951, 0.402464, 1.0)),
                    (1.0, (0.111678, 0.077768, 0.360944, 1.0)),
                ),
            )
        with g.Frame("coolwarm"):
            _string_46 = g.String(
                string="coolwarm: Diverging. 14 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_46 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.042129, 0.070875, 0.527059, 1.0)),
                    (0.078944, (0.084121, 0.155124, 0.723666, 1.0)),
                    (0.161856, (0.153989, 0.279601, 0.894763, 1.0)),
                    (0.244637, (0.256246, 0.422862, 0.991543, 1.0)),
                    (0.320233, (0.37837, 0.551122, 0.997868, 1.0)),
                    (0.405817, (0.53771, 0.667583, 0.911153, 1.0)),
                    (0.497386, (0.719273, 0.725601, 0.729085, 1.0)),
                    (0.58414, (0.871737, 0.623422, 0.50647, 1.0)),
                    (0.673169, (0.938203, 0.467399, 0.316952, 1.0)),
                    (0.762651, (0.901555, 0.297458, 0.177697, 1.0)),
                    (0.845734, (0.785211, 0.156212, 0.092642, 1.0)),
                    (0.922891, (0.628886, 0.059405, 0.045091, 1.0)),
                    (0.985277, (0.488686, 0.006856, 0.022971, 1.0)),
                    (1.0, (0.455946, 0.00111, 0.019925, 1.0)),
                ),
            )
        with g.Frame("bwr"):
            _string_47 = g.String(
                string="bwr: Diverging. 17 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_47 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 1.0, 1.0)),
                    (0.037333, (0.006115, 0.006115, 1.0, 1.0)),
                    (0.073352, (0.018106, 0.018106, 1.0, 1.0)),
                    (0.120206, (0.045677, 0.045677, 1.0, 1.0)),
                    (0.179364, (0.103175, 0.103175, 1.0, 1.0)),
                    (0.252703, (0.214956, 0.214956, 1.0, 1.0)),
                    (0.334631, (0.400159, 0.400156, 1.0, 1.0)),
                    (0.419759, (0.667293, 0.667322, 1.0, 1.0)),
                    (0.5, (0.995326, 0.994436, 0.995326, 1.0)),
                    (0.580241, (1.0, 0.667322, 0.667293, 1.0)),
                    (0.665369, (1.0, 0.400156, 0.400159, 1.0)),
                    (0.747297, (1.0, 0.214956, 0.214956, 1.0)),
                    (0.820636, (1.0, 0.103175, 0.103175, 1.0)),
                    (0.879794, (1.0, 0.045677, 0.045677, 1.0)),
                    (0.926648, (1.0, 0.018106, 0.018106, 1.0)),
                    (0.962667, (1.0, 0.006115, 0.006115, 1.0)),
                    (1.0, (1.0, 0.0, 0.0, 1.0)),
                ),
            )
        with g.Frame("seismic"):
            _string_48 = g.String(
                string="seismic: Diverging. 26 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_48 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.072208, 1.0)),
                    (0.038257, (0.0, 0.0, 0.135449, 1.0)),
                    (0.088056, (0.0, 0.0, 0.255832, 1.0)),
                    (0.141985, (0.0, 0.0, 0.440232, 1.0)),
                    (0.193369, (0.0, 0.0, 0.672584, 1.0)),
                    (0.251474, (0.000125, 0.000125, 1.0, 1.0)),
                    (0.271899, (0.007645, 0.007645, 0.999997, 1.0)),
                    (0.289177, (0.020425, 0.020425, 0.999998, 1.0)),
                    (0.311496, (0.04797, 0.04797, 0.999998, 1.0)),
                    (0.336586, (0.096261, 0.096261, 0.999999, 1.0)),
                    (0.371129, (0.195839, 0.195839, 1.0, 1.0)),
                    (0.413173, (0.377827, 0.377827, 1.0, 1.0)),
                    (0.456421, (0.642061, 0.642061, 1.0, 1.0)),
                    (0.498239, (0.978776, 0.978776, 1.0, 1.0)),
                    (0.501903, (1.0, 0.977631, 0.977631, 1.0)),
                    (0.542974, (1.0, 0.646936, 0.646936, 1.0)),
                    (0.58278, (1.0, 0.400142, 0.400142, 1.0)),
                    (0.621488, (1.0, 0.223368, 0.223368, 1.0)),
                    (0.656154, (1.0, 0.1138, 0.1138, 1.0)),
                    (0.684712, (1.0, 0.053671, 0.053671, 1.0)),
                    (0.710381, (1.0, 0.020744, 0.020744, 1.0)),
                    (0.728039, (1.0, 0.007695, 0.007695, 1.0)),
                    (0.748315, (1.0, 0.000145, 0.000145, 1.0)),
                    (0.834911, (0.650475, 0.0, 0.0, 1.0)),
                    (0.920288, (0.387094, 0.0, 0.0, 1.0)),
                    (1.0, (0.21077, 0.0, 0.0, 1.0)),
                ),
            )
        with g.Frame("berlin"):
            _string_49 = g.String(
                string="berlin: Diverging. 17 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_49 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.353122, 0.43412, 1.0, 1.0)),
                    (0.024332, (0.278108, 0.428053, 0.950624, 1.0)),
                    (0.064149, (0.189697, 0.410094, 0.850892, 1.0)),
                    (0.122923, (0.071917, 0.368818, 0.677606, 1.0)),
                    (0.19097, (0.03038, 0.22614, 0.394384, 1.0)),
                    (0.272548, (0.01709, 0.101879, 0.175422, 1.0)),
                    (0.36975, (0.005862, 0.023328, 0.03567, 1.0)),
                    (0.4678, (0.003671, 0.002372, 0.002385, 1.0)),
                    (0.553107, (0.01666, 0.003738, 0.0, 1.0)),
                    (0.641088, (0.054278, 0.006356, 0.0, 1.0)),
                    (0.715252, (0.127018, 0.012564, 0.002332, 1.0)),
                    (0.783759, (0.268713, 0.050139, 0.019554, 1.0)),
                    (0.850264, (0.422723, 0.115674, 0.079073, 1.0)),
                    (0.907499, (0.614019, 0.206657, 0.173389, 1.0)),
                    (0.951778, (0.772327, 0.29285, 0.270087, 1.0)),
                    (0.9874, (0.956935, 0.392214, 0.38627, 1.0)),
                    (1.0, (1.0, 0.422876, 0.423397, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("managua"):
            _string_50 = g.String(
                string="managua: Diverging. 14 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_50 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (1.0, 0.633739, 0.136724, 1.0)),
                    (0.020529, (0.956538, 0.559668, 0.127423, 1.0)),
                    (0.064285, (0.824503, 0.425558, 0.106838, 1.0)),
                    (0.126453, (0.684426, 0.283568, 0.08691, 1.0)),
                    (0.205888, (0.523251, 0.163196, 0.063131, 1.0)),
                    (0.309407, (0.330762, 0.066015, 0.046331, 1.0)),
                    (0.426404, (0.129498, 0.018523, 0.032772, 1.0)),
                    (0.556076, (0.0643, 0.019459, 0.080164, 1.0)),
                    (0.675879, (0.07466, 0.087741, 0.31291, 1.0)),
                    (0.782917, (0.10543, 0.219829, 0.521191, 1.0)),
                    (0.86407, (0.141771, 0.368716, 0.676177, 1.0)),
                    (0.928931, (0.174973, 0.546843, 0.815475, 1.0)),
                    (0.97573, (0.203784, 0.709257, 0.949283, 1.0)),
                    (1.0, (0.21962, 0.813587, 1.0, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("vanimo"):
            _string_51 = g.String(
                string="vanimo: Diverging. 17 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_51 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (1.0, 0.619979, 0.99213, 1.0)),
                    (0.016547, (0.957121, 0.533086, 0.908091, 1.0)),
                    (0.048742, (0.832104, 0.397504, 0.768565, 1.0)),
                    (0.094884, (0.7038, 0.248669, 0.598632, 1.0)),
                    (0.153217, (0.531133, 0.131501, 0.425183, 1.0)),
                    (0.222605, (0.371006, 0.062268, 0.277101, 1.0)),
                    (0.306888, (0.133746, 0.019284, 0.097896, 1.0)),
                    (0.394522, (0.018785, 0.007951, 0.014542, 1.0)),
                    (0.502054, (0.007271, 0.004775, 0.004662, 1.0)),
                    (0.614225, (0.012969, 0.018238, 0.006005, 1.0)),
                    (0.711858, (0.055506, 0.106691, 0.014496, 1.0)),
                    (0.795444, (0.115754, 0.231069, 0.025129, 1.0)),
                    (0.857759, (0.177102, 0.354955, 0.044844, 1.0)),
                    (0.909697, (0.266438, 0.526017, 0.086829, 1.0)),
                    (0.951054, (0.36968, 0.71394, 0.180724, 1.0)),
                    (0.979446, (0.445517, 0.852388, 0.273018, 1.0)),
                    (1.0, (0.524616, 1.0, 0.392334, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("twilight"):
            _string_52 = g.String(
                string="twilight: Cyclic. 22 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_52 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.775403, 0.69848, 0.777396, 1.0)),
                    (0.040866, (0.663365, 0.682931, 0.708375, 1.0)),
                    (0.096786, (0.371344, 0.531788, 0.583145, 1.0)),
                    (0.159571, (0.195103, 0.371151, 0.541975, 1.0)),
                    (0.220845, (0.128132, 0.231839, 0.520256, 1.0)),
                    (0.278209, (0.114394, 0.126217, 0.468463, 1.0)),
                    (0.33293, (0.110521, 0.052232, 0.388552, 1.0)),
                    (0.384314, (0.103611, 0.014139, 0.25235, 1.0)),
                    (0.432994, (0.060714, 0.00528, 0.103451, 1.0)),
                    (0.481784, (0.030346, 0.005766, 0.03891, 1.0)),
                    (0.507226, (0.027918, 0.007109, 0.037649, 1.0)),
                    (0.543128, (0.046717, 0.004446, 0.042477, 1.0)),
                    (0.593042, (0.101925, 0.008918, 0.071103, 1.0)),
                    (0.653165, (0.236826, 0.014992, 0.085048, 1.0)),
                    (0.713694, (0.377234, 0.048966, 0.075764, 1.0)),
                    (0.76835, (0.488077, 0.11785, 0.086922, 1.0)),
                    (0.808026, (0.539096, 0.194325, 0.118887, 1.0)),
                    (0.839906, (0.572367, 0.261296, 0.154623, 1.0)),
                    (0.873752, (0.603894, 0.370709, 0.247675, 1.0)),
                    (0.908279, (0.634832, 0.460435, 0.350861, 1.0)),
                    (0.956969, (0.729385, 0.648824, 0.633526, 1.0)),
                    (1.0, (0.767099, 0.704928, 0.788788, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("twilight_shifted"):
            _string_53 = g.String(
                string="twilight_shifted: Cyclic. 16 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_53 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.023373, 0.007044, 0.030903, 1.0)),
                    (0.061843, (0.05637, 0.003628, 0.080073, 1.0)),
                    (0.132771, (0.116136, 0.01734, 0.312704, 1.0)),
                    (0.206149, (0.109736, 0.099287, 0.460214, 1.0)),
                    (0.277526, (0.123859, 0.229365, 0.520702, 1.0)),
                    (0.343664, (0.204544, 0.381392, 0.53889, 1.0)),
                    (0.402378, (0.372262, 0.533747, 0.594638, 1.0)),
                    (0.458192, (0.658978, 0.671931, 0.690992, 1.0)),
                    (0.508873, (0.813745, 0.737213, 0.836921, 1.0)),
                    (0.58217, (0.628388, 0.502799, 0.376066, 1.0)),
                    (0.652464, (0.590523, 0.275941, 0.152815, 1.0)),
                    (0.727941, (0.492333, 0.110208, 0.078509, 1.0)),
                    (0.81303, (0.320152, 0.022988, 0.080926, 1.0)),
                    (0.894086, (0.115268, 0.008639, 0.080121, 1.0)),
                    (0.957345, (0.04118, 0.00463, 0.03959, 1.0)),
                    (1.0, (0.026915, 0.007218, 0.037617, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("hsv"):
            _string_54 = g.String(
                string="hsv: Cyclic. 32 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 1.7/255 in sRGB."
            )
            color_ramp_54 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (1.0, 0.0, 0.0, 1.0)),
                    (0.009057, (1.0, 0.002446, 0.0, 1.0)),
                    (0.035556, (1.0, 0.02216, 0.0, 1.0)),
                    (0.075376, (1.0, 0.137315, 0.0, 1.0)),
                    (0.119613, (0.995586, 0.441931, 0.0, 1.0)),
                    (0.151699, (1.0, 0.788854, 0.0, 1.0)),
                    (0.170053, (1.0, 1.0, 0.0, 1.0)),
                    (0.196125, (0.680775, 1.0, 0.0, 1.0)),
                    (0.236938, (0.28515, 1.0, 0.0, 1.0)),
                    (0.28335, (0.062993, 1.0, 0.0, 1.0)),
                    (0.321613, (0.0, 1.0, 0.0, 1.0)),
                    (0.353409, (0.000418, 1.0, 0.0, 1.0)),
                    (0.389132, (0.0, 1.0, 0.051941, 1.0)),
                    (0.42937, (0.0, 1.0, 0.234783, 1.0)),
                    (0.465892, (0.0, 1.0, 0.494528, 1.0)),
                    (0.497218, (0.0, 1.0, 0.893232, 1.0)),
                    (0.508543, (0.0, 1.0, 1.0, 1.0)),
                    (0.522777, (0.0, 0.834417, 1.0, 1.0)),
                    (0.557959, (0.0, 0.424571, 1.0, 1.0)),
                    (0.601444, (0.0, 0.146872, 1.0, 1.0)),
                    (0.643274, (0.0, 0.015269, 1.0, 1.0)),
                    (0.675379, (0.000292, 0.0, 1.0, 1.0)),
                    (0.686732, (0.003771, 0.0, 1.0, 1.0)),
                    (0.704516, (0.014526, 0.00015, 1.0, 1.0)),
                    (0.738668, (0.083059, 0.0, 1.0, 1.0)),
                    (0.785178, (0.326212, 0.0, 1.0, 1.0)),
                    (0.825701, (0.758258, 0.0, 1.0, 1.0)),
                    (0.847772, (1.0, 0.0, 1.0, 1.0)),
                    (0.875018, (1.0, 0.0, 0.657683, 1.0)),
                    (0.91664, (1.0, 0.0, 0.269041, 1.0)),
                    (0.962428, (1.0, 0.0, 0.058055, 1.0)),
                    (1.0, (1.0, 0.0, 0.0, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("Pastel1"):
            _string_55 = g.String(
                string="Pastel1: Qualitative. 9 colours copied exactly, constant interpolation."
            )
            color_ramp_55 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.964686, 0.456411, 0.423268, 1.0)),
                    (0.111111, (0.450786, 0.610496, 0.768151, 1.0)),
                    (0.222222, (0.603827, 0.83077, 0.55834, 1.0)),
                    (0.333333, (0.730461, 0.597202, 0.775822, 1.0)),
                    (0.444444, (0.991102, 0.693872, 0.381326, 1.0)),
                    (0.555556, (1.0, 1.0, 0.603827, 1.0)),
                    (0.666667, (0.783538, 0.686685, 0.508881, 1.0)),
                    (0.777778, (0.982251, 0.701102, 0.838799, 1.0)),
                    (0.888889, (0.887923, 0.887923, 0.887923, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("Pastel2"):
            _string_56 = g.String(
                string="Pastel2: Qualitative. 8 colours copied exactly, constant interpolation."
            )
            color_ramp_56 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.450786, 0.760525, 0.610496, 1.0)),
                    (0.125, (0.982251, 0.610496, 0.412543, 1.0)),
                    (0.25, (0.597202, 0.665387, 0.806952, 1.0)),
                    (0.375, (0.904661, 0.590619, 0.775822, 1.0)),
                    (0.5, (0.791298, 0.913099, 0.584078, 1.0)),
                    (0.625, (1.0, 0.887923, 0.423268, 1.0)),
                    (0.75, (0.879622, 0.760525, 0.603827, 1.0)),
                    (0.875, (0.603827, 0.603827, 0.603827, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("Paired"):
            _string_57 = g.String(
                string="Paired: Qualitative. 12 colours copied exactly, constant interpolation."
            )
            color_ramp_57 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.381326, 0.617207, 0.768151, 1.0)),
                    (0.083333, (0.013702, 0.187821, 0.456411, 1.0)),
                    (0.166667, (0.445201, 0.73791, 0.254152, 1.0)),
                    (0.25, (0.033105, 0.351533, 0.025187, 1.0)),
                    (0.333333, (0.964686, 0.323143, 0.318547, 1.0)),
                    (0.416667, (0.768151, 0.01033, 0.011612, 1.0)),
                    (0.5, (0.982251, 0.520996, 0.158961, 1.0)),
                    (0.583333, (1.0, 0.212231, 0.0, 1.0)),
                    (0.666667, (0.590619, 0.445201, 0.672443, 1.0)),
                    (0.75, (0.144128, 0.046665, 0.323143, 1.0)),
                    (0.833333, (1.0, 1.0, 0.318547, 1.0)),
                    (0.916667, (0.439657, 0.099899, 0.021219, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("Accent"):
            _string_58 = g.String(
                string="Accent: Qualitative. 8 colours copied exactly, constant interpolation."
            )
            color_ramp_58 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.212231, 0.584078, 0.212231, 1.0)),
                    (0.125, (0.514918, 0.423268, 0.658375, 1.0)),
                    (0.25, (0.982251, 0.527115, 0.238398, 1.0)),
                    (0.375, (1.0, 1.0, 0.318547, 1.0)),
                    (0.5, (0.039546, 0.14996, 0.434154, 1.0)),
                    (0.625, (0.871367, 0.000607, 0.212231, 1.0)),
                    (0.75, (0.520996, 0.104616, 0.008568, 1.0)),
                    (0.875, (0.132868, 0.132868, 0.132868, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("Dark2"):
            _string_59 = g.String(
                string="Dark2: Qualitative. 8 colours copied exactly, constant interpolation."
            )
            color_ramp_59 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.01096, 0.341914, 0.184475, 1.0)),
                    (0.125, (0.693872, 0.114435, 0.000607, 1.0)),
                    (0.25, (0.177888, 0.162029, 0.450786, 1.0)),
                    (0.375, (0.799103, 0.022174, 0.254152, 1.0)),
                    (0.5, (0.132868, 0.381326, 0.012983, 1.0)),
                    (0.625, (0.791298, 0.40724, 0.000607, 1.0)),
                    (0.75, (0.381326, 0.181164, 0.012286, 1.0)),
                    (0.875, (0.132868, 0.132868, 0.132868, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("Set1"):
            _string_60 = g.String(
                string="Set1: Qualitative. 9 colours copied exactly, constant interpolation."
            )
            color_ramp_60 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.775822, 0.01033, 0.011612, 1.0)),
                    (0.111111, (0.038204, 0.208637, 0.47932, 1.0)),
                    (0.222222, (0.074214, 0.42869, 0.068478, 1.0)),
                    (0.333333, (0.313989, 0.076185, 0.366253, 1.0)),
                    (0.444444, (1.0, 0.212231, 0.0, 1.0)),
                    (0.555556, (1.0, 1.0, 0.033105, 1.0)),
                    (0.666667, (0.381326, 0.093059, 0.021219, 1.0)),
                    (0.777778, (0.930111, 0.219526, 0.520996, 1.0)),
                    (0.888889, (0.318547, 0.318547, 0.318547, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("Set2"):
            _string_61 = g.String(
                string="Set2: Qualitative. 8 colours copied exactly, constant interpolation."
            )
            color_ramp_61 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.132868, 0.539479, 0.376262, 1.0)),
                    (0.125, (0.973445, 0.266356, 0.122139, 1.0)),
                    (0.25, (0.266356, 0.351533, 0.597202, 1.0)),
                    (0.375, (0.799103, 0.254152, 0.545724, 1.0)),
                    (0.5, (0.381326, 0.686685, 0.088656, 1.0)),
                    (0.625, (1.0, 0.693872, 0.028426, 1.0)),
                    (0.75, (0.783538, 0.552011, 0.296138, 1.0)),
                    (0.875, (0.450786, 0.450786, 0.450786, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("Set3"):
            _string_62 = g.String(
                string="Set3: Qualitative. 12 colours copied exactly, constant interpolation."
            )
            color_ramp_62 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.266356, 0.651406, 0.571125, 1.0)),
                    (0.083333, (1.0, 1.0, 0.450786, 1.0)),
                    (0.166667, (0.514918, 0.491021, 0.701102, 1.0)),
                    (0.25, (0.964686, 0.215861, 0.168269, 1.0)),
                    (0.333333, (0.215861, 0.439657, 0.651406, 1.0)),
                    (0.416667, (0.982251, 0.456411, 0.122139, 1.0)),
                    (0.5, (0.450786, 0.730461, 0.141263, 1.0)),
                    (0.583333, (0.973445, 0.610496, 0.783538, 1.0)),
                    (0.666667, (0.693872, 0.693872, 0.693872, 1.0)),
                    (0.75, (0.502886, 0.215861, 0.508881, 1.0)),
                    (0.833333, (0.603827, 0.83077, 0.55834, 1.0)),
                    (0.916667, (1.0, 0.846873, 0.158961, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("tab10"):
            _string_63 = g.String(
                string="tab10: Qualitative. 10 colours copied exactly, constant interpolation."
            )
            color_ramp_63 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.013702, 0.184475, 0.456411, 1.0)),
                    (0.1, (1.0, 0.212231, 0.004391, 1.0)),
                    (0.2, (0.025187, 0.351533, 0.025187, 1.0)),
                    (0.3, (0.672443, 0.020289, 0.021219, 1.0)),
                    (0.4, (0.296138, 0.135633, 0.508881, 1.0)),
                    (0.5, (0.262251, 0.093059, 0.07036, 1.0)),
                    (0.6, (0.768151, 0.184475, 0.539479, 1.0)),
                    (0.7, (0.212231, 0.212231, 0.212231, 1.0)),
                    (0.8, (0.502886, 0.508881, 0.015996, 1.0)),
                    (0.9, (0.008568, 0.514918, 0.62396, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("tab20"):
            _string_64 = g.String(
                string="tab20: Qualitative. 20 colours copied exactly, constant interpolation."
            )
            color_ramp_64 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.013702, 0.184475, 0.456411, 1.0)),
                    (0.05, (0.423268, 0.571125, 0.806952, 1.0)),
                    (0.1, (1.0, 0.212231, 0.004391, 1.0)),
                    (0.15, (1.0, 0.496933, 0.187821, 1.0)),
                    (0.2, (0.025187, 0.351533, 0.025187, 1.0)),
                    (0.25, (0.313989, 0.73791, 0.254152, 1.0)),
                    (0.3, (0.672443, 0.020289, 0.021219, 1.0)),
                    (0.35, (1.0, 0.313989, 0.304987, 1.0)),
                    (0.4, (0.296138, 0.135633, 0.508881, 1.0)),
                    (0.45, (0.55834, 0.434154, 0.665387, 1.0)),
                    (0.5, (0.262251, 0.093059, 0.07036, 1.0)),
                    (0.55, (0.552011, 0.332452, 0.296138, 1.0)),
                    (0.6, (0.768151, 0.184475, 0.539479, 1.0)),
                    (0.65, (0.930111, 0.467784, 0.64448, 1.0)),
                    (0.7, (0.212231, 0.212231, 0.212231, 1.0)),
                    (0.75, (0.571125, 0.571125, 0.571125, 1.0)),
                    (0.8, (0.502886, 0.508881, 0.015996, 1.0)),
                    (0.85, (0.708376, 0.708376, 0.266356, 1.0)),
                    (0.9, (0.008568, 0.514918, 0.62396, 1.0)),
                    (0.95, (0.341914, 0.701102, 0.783538, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("tab20b"):
            _string_65 = g.String(
                string="tab20b: Qualitative. 20 colours copied exactly, constant interpolation."
            )
            color_ramp_65 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.040915, 0.043735, 0.191202, 1.0)),
                    (0.05, (0.084376, 0.088656, 0.366253, 1.0)),
                    (0.1, (0.147027, 0.155926, 0.62396, 1.0)),
                    (0.15, (0.332452, 0.341914, 0.730461, 1.0)),
                    (0.2, (0.124772, 0.191202, 0.040915, 1.0)),
                    (0.25, (0.262251, 0.361307, 0.084376, 1.0)),
                    (0.3, (0.462077, 0.62396, 0.147027, 1.0)),
                    (0.35, (0.617207, 0.708376, 0.332452, 1.0)),
                    (0.4, (0.262251, 0.152926, 0.030713, 1.0)),
                    (0.45, (0.508881, 0.341914, 0.040915, 1.0)),
                    (0.5, (0.799103, 0.491021, 0.084376, 1.0)),
                    (0.55, (0.799103, 0.597202, 0.296138, 1.0)),
                    (0.6, (0.23074, 0.045186, 0.040915, 1.0)),
                    (0.65, (0.417885, 0.066626, 0.068478, 1.0)),
                    (0.7, (0.672443, 0.119538, 0.147027, 1.0)),
                    (0.75, (0.799103, 0.304987, 0.332452, 1.0)),
                    (0.8, (0.198069, 0.052861, 0.171441, 1.0)),
                    (0.85, (0.376262, 0.082283, 0.296138, 1.0)),
                    (0.9, (0.617207, 0.152926, 0.508881, 1.0)),
                    (0.95, (0.730461, 0.341914, 0.672443, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("tab20c"):
            _string_66 = g.String(
                string="tab20c: Qualitative. 20 colours copied exactly, constant interpolation."
            )
            color_ramp_66 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.030713, 0.223228, 0.508881, 1.0)),
                    (0.05, (0.147027, 0.423268, 0.672443, 1.0)),
                    (0.1, (0.341914, 0.590619, 0.752942, 1.0)),
                    (0.15, (0.564712, 0.708376, 0.863157, 1.0)),
                    (0.2, (0.791298, 0.090842, 0.004025, 1.0)),
                    (0.25, (0.982251, 0.266356, 0.045186, 1.0)),
                    (0.3, (0.982251, 0.423268, 0.147027, 1.0)),
                    (0.35, (0.982251, 0.630757, 0.361307, 1.0)),
                    (0.4, (0.030713, 0.366253, 0.088656, 1.0)),
                    (0.45, (0.174647, 0.552011, 0.181164, 1.0)),
                    (0.5, (0.3564, 0.693872, 0.327778, 1.0)),
                    (0.55, (0.571125, 0.814847, 0.527115, 1.0)),
                    (0.6, (0.177888, 0.147027, 0.439657, 1.0)),
                    (0.65, (0.341914, 0.323143, 0.57758, 1.0)),
                    (0.7, (0.502886, 0.508881, 0.715694, 1.0)),
                    (0.75, (0.701102, 0.701102, 0.83077, 1.0)),
                    (0.8, (0.124772, 0.124772, 0.124772, 1.0)),
                    (0.85, (0.304987, 0.304987, 0.304987, 1.0)),
                    (0.9, (0.508881, 0.508881, 0.508881, 1.0)),
                    (0.95, (0.693872, 0.693872, 0.693872, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("okabe_ito"):
            _string_67 = g.String(
                string="okabe_ito: Qualitative. 8 colours copied exactly, constant interpolation."
            )
            color_ramp_67 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.0, 1.0)),
                    (0.125, (0.791298, 0.346704, 0.0, 1.0)),
                    (0.25, (0.093059, 0.456411, 0.814847, 1.0)),
                    (0.375, (0.0, 0.341914, 0.171441, 1.0)),
                    (0.5, (0.871367, 0.775822, 0.05448, 1.0)),
                    (0.625, (0.0, 0.168269, 0.445201, 1.0)),
                    (0.75, (0.665387, 0.111932, 0.0, 1.0)),
                    (0.875, (0.603827, 0.191202, 0.386429, 1.0)),
                ),
                color_interpolation="CONSTANT",
            )
        with g.Frame("ocean"):
            _string_68 = g.String(
                string="ocean: Miscellaneous. 19 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_68 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.211977, 0.0, 1.0)),
                    (0.077069, (0.0, 0.120016, 0.006371, 1.0)),
                    (0.155622, (0.0, 0.055917, 0.020118, 1.0)),
                    (0.228374, (0.0, 0.020324, 0.042194, 1.0)),
                    (0.28298, (0.0, 0.006221, 0.064729, 1.0)),
                    (0.33257, (0.0, 0.0, 0.090076, 1.0)),
                    (0.380001, (0.0, 0.005596, 0.118958, 1.0)),
                    (0.433868, (0.0, 0.018815, 0.157196, 1.0)),
                    (0.503393, (0.0, 0.051272, 0.216277, 1.0)),
                    (0.582013, (0.0, 0.112437, 0.296746, 1.0)),
                    (0.66729, (0.0, 0.212127, 0.401493, 1.0)),
                    (0.692389, (0.006397, 0.253177, 0.437977, 1.0)),
                    (0.716807, (0.018827, 0.289001, 0.47169, 1.0)),
                    (0.750387, (0.049753, 0.349602, 0.523212, 1.0)),
                    (0.787307, (0.105458, 0.420379, 0.58213, 1.0)),
                    (0.834146, (0.212117, 0.523139, 0.662725, 1.0)),
                    (0.890824, (0.404184, 0.665639, 0.76873, 1.0)),
                    (0.946484, (0.667271, 0.825558, 0.881897, 1.0)),
                    (1.0, (0.995082, 0.998805, 0.999474, 1.0)),
                ),
            )
        with g.Frame("gist_earth"):
            _string_69 = g.String(
                string="gist_earth: Miscellaneous. 23 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_69 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.0, 1.0)),
                    (0.003618, (0.000186, 0.0, 0.023038, 1.0)),
                    (0.014828, (0.000765, 0.0, 0.068885, 1.0)),
                    (0.027462, (0.001416, 0.000014, 0.154484, 1.0)),
                    (0.031424, (0.001626, 0.000683, 0.174682, 1.0)),
                    (0.057443, (0.002954, 0.005526, 0.178405, 1.0)),
                    (0.086113, (0.004571, 0.015685, 0.183392, 1.0)),
                    (0.134923, (0.008393, 0.046222, 0.189702, 1.0)),
                    (0.19942, (0.015634, 0.110149, 0.201918, 1.0)),
                    (0.282353, (0.029163, 0.213386, 0.212477, 1.0)),
                    (0.373595, (0.043595, 0.269361, 0.123353, 1.0)),
                    (0.457532, (0.059369, 0.320377, 0.065241, 1.0)),
                    (0.469437, (0.069457, 0.330626, 0.06114, 1.0)),
                    (0.511473, (0.124437, 0.359968, 0.075075, 1.0)),
                    (0.623226, (0.309568, 0.426336, 0.09899, 1.0)),
                    (0.698875, (0.472458, 0.472539, 0.111213, 1.0)),
                    (0.786223, (0.528633, 0.365712, 0.131707, 1.0)),
                    (0.820598, (0.588422, 0.389124, 0.198875, 1.0)),
                    (0.86468, (0.675521, 0.444809, 0.308164, 1.0)),
                    (0.897402, (0.743034, 0.508418, 0.410873, 1.0)),
                    (0.937981, (0.832897, 0.638817, 0.596378, 1.0)),
                    (0.966142, (0.898942, 0.757018, 0.75082, 1.0)),
                    (1.0, (0.982009, 0.963496, 0.962591, 1.0)),
                ),
            )
        with g.Frame("terrain"):
            _string_70 = g.String(
                string="terrain: Miscellaneous. 26 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_70 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.032983, 0.032854, 0.318175, 1.0)),
                    (0.022427, (0.024518, 0.054136, 0.39179, 1.0)),
                    (0.054151, (0.014524, 0.095738, 0.511937, 1.0)),
                    (0.097053, (0.00577, 0.175107, 0.704689, 1.0)),
                    (0.149499, (0.0, 0.313506, 0.992992, 1.0)),
                    (0.185473, (0.0, 0.409647, 0.575027, 1.0)),
                    (0.213381, (0.0, 0.484972, 0.337448, 1.0)),
                    (0.236496, (0.0, 0.560207, 0.194145, 1.0)),
                    (0.25026, (0.0, 0.603438, 0.132092, 1.0)),
                    (0.269162, (0.006273, 0.630629, 0.1445, 1.0)),
                    (0.287991, (0.019244, 0.656108, 0.154673, 1.0)),
                    (0.309156, (0.044318, 0.687843, 0.16875, 1.0)),
                    (0.339277, (0.10198, 0.730857, 0.188246, 1.0)),
                    (0.375066, (0.21009, 0.789115, 0.214019, 1.0)),
                    (0.417109, (0.398573, 0.852891, 0.24628, 1.0)),
                    (0.460777, (0.67282, 0.934276, 0.283124, 1.0)),
                    (0.499509, (0.994766, 0.992088, 0.317593, 1.0)),
                    (0.581748, (0.663074, 0.58072, 0.223905, 1.0)),
                    (0.651624, (0.440357, 0.328008, 0.158787, 1.0)),
                    (0.706913, (0.300926, 0.184982, 0.117282, 1.0)),
                    (0.750169, (0.212909, 0.104916, 0.087742, 1.0)),
                    (0.79114, (0.296766, 0.181134, 0.160465, 1.0)),
                    (0.841075, (0.420864, 0.307291, 0.285535, 1.0)),
                    (0.895497, (0.586257, 0.49176, 0.472674, 1.0)),
                    (0.949972, (0.78493, 0.729138, 0.717432, 1.0)),
                    (1.0, (0.998189, 0.997015, 0.996724, 1.0)),
                ),
            )
        with g.Frame("gist_stern"):
            _string_71 = g.String(
                string="gist_stern: Miscellaneous. 32 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.8/255 in sRGB."
            )
            color_ramp_71 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.0, 1.0)),
                    (0.003081, (0.003432, 0.000246, 0.000528, 1.0)),
                    (0.009975, (0.025072, 0.000761, 0.001463, 1.0)),
                    (0.017805, (0.082901, 0.001406, 0.002858, 1.0)),
                    (0.02612, (0.186414, 0.001981, 0.003935, 1.0)),
                    (0.038112, (0.427403, 0.00301, 0.006946, 1.0)),
                    (0.054902, (0.982356, 0.00423, 0.010903, 1.0)),
                    (0.107164, (0.49383, 0.01084, 0.036062, 1.0)),
                    (0.146045, (0.253684, 0.018775, 0.069825, 1.0)),
                    (0.173659, (0.135661, 0.025324, 0.097993, 1.0)),
                    (0.193795, (0.075058, 0.031283, 0.124938, 1.0)),
                    (0.213811, (0.034203, 0.037447, 0.152253, 1.0)),
                    (0.230287, (0.013774, 0.043376, 0.17973, 1.0)),
                    (0.247038, (0.002112, 0.04966, 0.208306, 1.0)),
                    (0.247867, (0.048901, 0.048901, 0.204461, 1.0)),
                    (0.341672, (0.094319, 0.094319, 0.418016, 1.0)),
                    (0.414463, (0.142329, 0.142329, 0.651126, 1.0)),
                    (0.499156, (0.212002, 0.212002, 0.986875, 1.0)),
                    (0.565708, (0.27977, 0.27977, 0.465786, 1.0)),
                    (0.615423, (0.336358, 0.336358, 0.217302, 1.0)),
                    (0.652078, (0.382779, 0.382781, 0.099311, 1.0)),
                    (0.679681, (0.419403, 0.419383, 0.043636, 1.0)),
                    (0.703062, (0.452241, 0.452247, 0.015624, 1.0)),
                    (0.719571, (0.476196, 0.476196, 0.005077, 1.0)),
                    (0.735573, (0.500051, 0.500051, 0.0, 1.0)),
                    (0.755732, (0.53115, 0.53115, 0.006477, 1.0)),
                    (0.776937, (0.565116, 0.565116, 0.020598, 1.0)),
                    (0.80659, (0.614601, 0.614601, 0.057204, 1.0)),
                    (0.837954, (0.669688, 0.669688, 0.122009, 1.0)),
                    (0.878174, (0.744602, 0.744602, 0.248753, 1.0)),
                    (0.930995, (0.849217, 0.849217, 0.495368, 1.0)),
                    (1.0, (0.998979, 0.998979, 0.985029, 1.0)),
                ),
            )
        with g.Frame("gnuplot"):
            _string_72 = g.String(
                string="gnuplot: Miscellaneous. 30 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.7/255 in sRGB."
            )
            color_ramp_72 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.000009, 0.0, 0.0, 1.0)),
                    (0.00405, (0.005342, 0.0, 0.002024, 1.0)),
                    (0.010918, (0.010756, 0.000001, 0.005242, 1.0)),
                    (0.024646, (0.021219, 0.000001, 0.019577, 1.0)),
                    (0.037655, (0.031328, 0.000012, 0.043771, 1.0)),
                    (0.052771, (0.043025, 0.000001, 0.08515, 1.0)),
                    (0.069828, (0.056908, 0.000047, 0.148994, 1.0)),
                    (0.095259, (0.077298, 0.000034, 0.274762, 1.0)),
                    (0.182437, (0.152256, 0.000435, 0.817211, 1.0)),
                    (0.221556, (0.188181, 0.000869, 0.972615, 1.0)),
                    (0.25113, (0.214885, 0.001169, 1.0, 1.0)),
                    (0.280424, (0.242519, 0.001789, 0.967343, 1.0)),
                    (0.319222, (0.279036, 0.002377, 0.807949, 1.0)),
                    (0.403006, (0.360527, 0.00526, 0.284814, 1.0)),
                    (0.428244, (0.385902, 0.007152, 0.157382, 1.0)),
                    (0.445287, (0.402738, 0.008205, 0.091582, 1.0)),
                    (0.460486, (0.418074, 0.009526, 0.048055, 1.0)),
                    (0.473061, (0.430763, 0.011185, 0.023375, 1.0)),
                    (0.484662, (0.442511, 0.011939, 0.008792, 1.0)),
                    (0.496853, (0.455021, 0.014433, 0.001, 1.0)),
                    (0.503046, (0.461169, 0.014122, 0.0, 1.0)),
                    (0.577603, (0.538129, 0.029956, 0.0, 1.0)),
                    (0.641814, (0.605669, 0.055568, 0.0, 1.0)),
                    (0.703912, (0.671962, 0.097896, 0.0, 1.0)),
                    (0.763202, (0.736129, 0.16366, 0.0, 1.0)),
                    (0.820265, (0.798635, 0.262197, 0.0, 1.0)),
                    (0.867761, (0.851203, 0.381773, 0.0, 1.0)),
                    (0.908607, (0.896781, 0.519575, 0.0, 1.0)),
                    (0.95232, (0.945999, 0.712344, 0.0, 1.0)),
                    (1.0, (0.999944, 0.994345, 0.0, 1.0)),
                ),
            )
        with g.Frame("gnuplot2"):
            _string_73 = g.String(
                string="gnuplot2: Miscellaneous. 32 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 1.2/255 in sRGB."
            )
            color_ramp_73 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.0, 1.0)),
                    (0.008667, (0.0, 0.0, 0.001664, 1.0)),
                    (0.036881, (0.0, 0.0, 0.012573, 1.0)),
                    (0.07964, (0.0, 0.0, 0.068348, 1.0)),
                    (0.131947, (0.0, 0.0, 0.220751, 1.0)),
                    (0.179277, (0.0, 0.0, 0.475159, 1.0)),
                    (0.214636, (0.0, 0.0, 0.684013, 1.0)),
                    (0.245429, (0.0, 0.0, 0.995878, 1.0)),
                    (0.2562, (0.0, 0.0, 1.0, 1.0)),
                    (0.298876, (0.012278, 0.0, 1.0, 1.0)),
                    (0.355804, (0.077981, 0.0, 1.0, 1.0)),
                    (0.402661, (0.198815, 0.0, 1.0, 1.0)),
                    (0.428614, (0.260154, 0.0, 1.0, 1.0)),
                    (0.467685, (0.418725, 0.008031, 0.771771, 1.0)),
                    (0.511506, (0.60849, 0.023956, 0.648102, 1.0)),
                    (0.55572, (0.930709, 0.06333, 0.465256, 1.0)),
                    (0.570209, (1.0, 0.068689, 0.461886, 1.0)),
                    (0.612079, (1.0, 0.116837, 0.331635, 1.0)),
                    (0.670629, (1.0, 0.204734, 0.211553, 1.0)),
                    (0.746497, (1.0, 0.369596, 0.088235, 1.0)),
                    (0.822024, (1.0, 0.61485, 0.025532, 1.0)),
                    (0.876676, (1.0, 0.807857, 0.006415, 1.0)),
                    (0.915085, (1.0, 1.0, 0.0, 1.0)),
                    (0.921703, (1.0, 1.0, 0.0, 1.0)),
                    (0.924288, (1.0, 0.996958, 0.009142, 1.0)),
                    (0.928624, (1.0, 1.0, 0.000024, 1.0)),
                    (0.943922, (1.0, 1.0, 0.061631, 1.0)),
                    (0.960065, (1.0, 1.0, 0.1958, 1.0)),
                    (0.975327, (1.0, 1.0, 0.432812, 1.0)),
                    (0.987623, (1.0, 1.0, 0.665287, 1.0)),
                    (0.99823, (1.0, 1.0, 0.973958, 1.0)),
                    (1.0, (1.0, 1.0, 1.0, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("CMRmap"):
            _string_74 = g.String(
                string="CMRmap: Miscellaneous. 25 stops fitted to matplotlib 3.11.2, linear interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_74 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.0, 1.0)),
                    (0.018299, (0.001728, 0.001733, 0.005885, 1.0)),
                    (0.037996, (0.003487, 0.00348, 0.019005, 1.0)),
                    (0.064302, (0.006706, 0.006725, 0.052033, 1.0)),
                    (0.092673, (0.011812, 0.011713, 0.110882, 1.0)),
                    (0.125692, (0.019411, 0.01955, 0.213234, 1.0)),
                    (0.182863, (0.038789, 0.019545, 0.334482, 1.0)),
                    (0.250711, (0.072435, 0.019721, 0.520304, 1.0)),
                    (0.304535, (0.152212, 0.024773, 0.364718, 1.0)),
                    (0.371406, (0.303084, 0.032902, 0.218742, 1.0)),
                    (0.420684, (0.510837, 0.038534, 0.111991, 1.0)),
                    (0.461604, (0.741546, 0.045839, 0.052131, 1.0)),
                    (0.499931, (0.996298, 0.049356, 0.018703, 1.0)),
                    (0.561305, (0.893915, 0.111968, 0.006417, 1.0)),
                    (0.62268, (0.789052, 0.20726, 0.000011, 1.0)),
                    (0.687298, (0.786377, 0.344885, 0.00376, 1.0)),
                    (0.751083, (0.788111, 0.522661, 0.009836, 1.0)),
                    (0.778019, (0.787178, 0.577557, 0.028942, 1.0)),
                    (0.806091, (0.786859, 0.632597, 0.062199, 1.0)),
                    (0.837437, (0.788321, 0.700631, 0.117332, 1.0)),
                    (0.872979, (0.78604, 0.783007, 0.205529, 1.0)),
                    (0.897376, (0.822197, 0.823859, 0.303828, 1.0)),
                    (0.931068, (0.8797, 0.878566, 0.479296, 1.0)),
                    (0.966409, (0.939376, 0.940127, 0.716809, 1.0)),
                    (1.0, (0.999996, 0.999732, 0.996843, 1.0)),
                ),
            )
        with g.Frame("cubehelix"):
            _string_75 = g.String(
                string="cubehelix: Miscellaneous. 15 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.8/255 in sRGB."
            )
            color_ramp_75 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.000057, 0.0, 1.0)),
                    (0.01632, (0.00143, 0.000317, 0.001236, 1.0)),
                    (0.063822, (0.008078, 0.003117, 0.008847, 1.0)),
                    (0.128002, (0.012309, 0.008871, 0.043248, 1.0)),
                    (0.213756, (0.006512, 0.044524, 0.095524, 1.0)),
                    (0.315489, (0.007644, 0.154615, 0.041794, 1.0)),
                    (0.412835, (0.069517, 0.207976, 0.020457, 1.0)),
                    (0.508257, (0.373411, 0.187218, 0.047655, 1.0)),
                    (0.603913, (0.705038, 0.190455, 0.247716, 1.0)),
                    (0.707848, (0.630314, 0.325976, 0.794752, 1.0)),
                    (0.799843, (0.497082, 0.594575, 0.932399, 1.0)),
                    (0.878267, (0.585757, 0.822901, 0.862538, 1.0)),
                    (0.936794, (0.741124, 0.922142, 0.84503, 1.0)),
                    (0.983865, (0.957197, 0.986872, 0.960541, 1.0)),
                    (1.0, (1.0, 1.0, 1.0, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("brg"):
            _string_76 = g.String(
                string="brg: Miscellaneous. 20 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_76 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 1.0, 1.0)),
                    (0.016283, (0.001599, 0.0, 0.947795, 1.0)),
                    (0.063028, (0.010722, 0.0, 0.728073, 1.0)),
                    (0.123951, (0.043759, 0.0, 0.525747, 1.0)),
                    (0.201288, (0.121182, 0.0, 0.304456, 1.0)),
                    (0.290893, (0.288493, 0.0, 0.13305, 1.0)),
                    (0.370618, (0.502277, 0.0, 0.047066, 1.0)),
                    (0.43475, (0.733028, 0.0, 0.011465, 1.0)),
                    (0.476975, (0.903133, 0.0, 0.003044, 1.0)),
                    (0.496922, (1.0, 0.000142, 0.000169, 1.0)),
                    (0.510395, (0.982443, 0.000308, 0.0, 1.0)),
                    (0.557644, (0.738354, 0.010036, 0.0, 1.0)),
                    (0.608236, (0.582538, 0.034134, 0.0, 1.0)),
                    (0.672264, (0.375804, 0.088989, 0.0, 1.0)),
                    (0.743986, (0.220088, 0.19427, 0.0, 1.0)),
                    (0.821925, (0.094259, 0.360995, 0.0, 1.0)),
                    (0.892384, (0.032648, 0.580282, 0.0, 1.0)),
                    (0.947078, (0.008446, 0.770216, 0.0, 1.0)),
                    (0.987515, (0.001218, 0.961298, 0.0, 1.0)),
                    (1.0, (0.0, 1.0, 0.0, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("gist_rainbow"):
            _string_77 = g.String(
                string="gist_rainbow: Miscellaneous. 32 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 1.9/255 in sRGB."
            )
            color_ramp_77 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (1.0, 0.0, 0.024755, 1.0)),
                    (0.019361, (1.0, 0.0, 0.0, 1.0)),
                    (0.042276, (1.0, 0.0, 0.0, 1.0)),
                    (0.077039, (1.0, 0.036476, 0.0, 1.0)),
                    (0.119941, (1.0, 0.178697, 0.0, 1.0)),
                    (0.163874, (1.0, 0.467358, 0.0, 1.0)),
                    (0.196847, (1.0, 0.802627, 0.0, 1.0)),
                    (0.213277, (1.0, 1.0, 0.0, 1.0)),
                    (0.224658, (0.924999, 1.0, 0.0, 1.0)),
                    (0.263253, (0.467939, 1.0, 0.0, 1.0)),
                    (0.31021, (0.181342, 1.0, 0.0, 1.0)),
                    (0.355765, (0.0278, 1.0, 0.0, 1.0)),
                    (0.39096, (0.0, 1.0, 0.0, 1.0)),
                    (0.407916, (0.0, 1.0, 0.0, 1.0)),
                    (0.436573, (0.0, 1.0, 0.020787, 1.0)),
                    (0.479017, (0.0, 1.0, 0.121022, 1.0)),
                    (0.527404, (0.0, 1.0, 0.413868, 1.0)),
                    (0.56339, (0.0, 1.0, 0.743388, 1.0)),
                    (0.583485, (0.0, 1.0, 1.0, 1.0)),
                    (0.593282, (0.0, 0.948481, 1.0, 1.0)),
                    (0.62414, (0.0, 0.562378, 1.0, 1.0)),
                    (0.660049, (0.0, 0.30651, 1.0, 1.0)),
                    (0.701052, (0.0, 0.094395, 1.0, 1.0)),
                    (0.735063, (0.0, 0.023876, 1.0, 1.0)),
                    (0.759911, (0.0, 0.0, 1.0, 1.0)),
                    (0.783879, (0.0, 0.0, 1.0, 1.0)),
                    (0.823366, (0.044629, 0.0, 1.0, 1.0)),
                    (0.870305, (0.243474, 0.0, 1.0, 1.0)),
                    (0.912922, (0.52636, 0.0, 1.0, 1.0)),
                    (0.947887, (0.976989, 0.0, 1.0, 1.0)),
                    (0.960538, (1.0, 0.0, 1.0, 1.0)),
                    (1.0, (1.0, 0.0, 0.427911, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("rainbow"):
            _string_78 = g.String(
                string="rainbow: Miscellaneous. 24 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.0/255 in sRGB."
            )
            color_ramp_78 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.213828, 0.000033, 1.0, 1.0)),
                    (0.018301, (0.181594, 0.004459, 0.999233, 1.0)),
                    (0.034659, (0.154902, 0.010952, 0.996851, 1.0)),
                    (0.061647, (0.11683, 0.029578, 0.989777, 1.0)),
                    (0.09988, (0.07206, 0.075429, 0.973085, 1.0)),
                    (0.146728, (0.034105, 0.163742, 0.941675, 1.0)),
                    (0.189234, (0.013104, 0.273755, 0.903618, 1.0)),
                    (0.220274, (0.004613, 0.364063, 0.870639, 1.0)),
                    (0.250592, (0.0, 0.4603, 0.835041, 1.0)),
                    (0.285755, (0.005803, 0.572943, 0.789544, 1.0)),
                    (0.32457, (0.018447, 0.697464, 0.735093, 1.0)),
                    (0.377332, (0.051094, 0.843878, 0.655262, 1.0)),
                    (0.432946, (0.108285, 0.95605, 0.566509, 1.0)),
                    (0.488157, (0.1907, 1.0, 0.477326, 1.0)),
                    (0.546599, (0.307877, 0.98261, 0.384055, 1.0)),
                    (0.612018, (0.479761, 0.870909, 0.286528, 1.0)),
                    (0.67717, (0.69629, 0.693515, 0.199932, 1.0)),
                    (0.751249, (1.0, 0.451178, 0.118018, 1.0)),
                    (0.826008, (1.0, 0.229458, 0.057986, 1.0)),
                    (0.879718, (1.0, 0.10977, 0.028877, 1.0)),
                    (0.918679, (1.0, 0.050695, 0.014515, 1.0)),
                    (0.948327, (1.0, 0.021322, 0.00712, 1.0)),
                    (0.976818, (1.0, 0.005827, 0.002704, 1.0)),
                    (1.0, (1.0, 0.0, 0.000067, 1.0)),
                ),
            )
        with g.Frame("jet"):
            _string_79 = g.String(
                string="jet: Miscellaneous. 32 stops fitted to matplotlib 3.11.2, linear interpolation, max error 1.5/255 in sRGB."
            )
            color_ramp_79 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.21081, 1.0)),
                    (0.036254, (0.0, 0.0, 0.39379, 1.0)),
                    (0.072758, (0.0, 0.0, 0.651197, 1.0)),
                    (0.110597, (0.0, 0.0, 0.999979, 1.0)),
                    (0.125806, (0.0, 0.000084, 1.0, 1.0)),
                    (0.143901, (0.0, 0.006193, 1.0, 1.0)),
                    (0.163017, (0.0, 0.019144, 1.0, 1.0)),
                    (0.186991, (0.0, 0.04841, 1.0, 1.0)),
                    (0.217227, (0.0, 0.109549, 1.0, 1.0)),
                    (0.250187, (0.0, 0.211012, 1.0, 1.0)),
                    (0.292138, (0.0, 0.398271, 1.0, 1.0)),
                    (0.338938, (0.0, 0.695078, 1.0, 1.0)),
                    (0.350766, (0.0, 0.797221, 0.926333, 1.0)),
                    (0.375301, (0.006843, 0.996562, 0.75688, 1.0)),
                    (0.402452, (0.023071, 1.0, 0.598879, 1.0)),
                    (0.437665, (0.062752, 1.0, 0.424242, 1.0)),
                    (0.481029, (0.145562, 1.0, 0.25341, 1.0)),
                    (0.532265, (0.299451, 1.0, 0.114737, 1.0)),
                    (0.581269, (0.512108, 1.0, 0.037908, 1.0)),
                    (0.616664, (0.71001, 1.0, 0.010325, 1.0)),
                    (0.640226, (0.857055, 0.995777, 0.002133, 1.0)),
                    (0.651981, (0.948032, 0.907985, 0.0, 1.0)),
                    (0.661592, (1.0, 0.816026, 0.0, 1.0)),
                    (0.72341, (1.0, 0.426414, 0.0, 1.0)),
                    (0.777041, (1.0, 0.202197, 0.0, 1.0)),
                    (0.817083, (1.0, 0.094455, 0.0, 1.0)),
                    (0.848311, (1.0, 0.041294, 0.0, 1.0)),
                    (0.871743, (1.0, 0.017105, 0.0, 1.0)),
                    (0.891517, (0.99235, 0.005253, 0.0, 1.0)),
                    (0.911074, (0.785416, 0.0, 0.0, 1.0)),
                    (0.958761, (0.422013, 0.0, 0.0, 1.0)),
                    (1.0, (0.210086, 0.0, 0.0, 1.0)),
                ),
            )
        with g.Frame("turbo"):
            _string_80 = g.String(
                string="turbo: Miscellaneous. 20 stops fitted to matplotlib 3.11.2, b-spline interpolation, max error 0.9/255 in sRGB."
            )
            color_ramp_80 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.028687, 0.004685, 0.032998, 1.0)),
                    (0.017026, (0.037297, 0.013856, 0.10371, 1.0)),
                    (0.037051, (0.041999, 0.024247, 0.185562, 1.0)),
                    (0.068966, (0.054588, 0.055981, 0.393414, 1.0)),
                    (0.110247, (0.060452, 0.112387, 0.67523, 1.0)),
                    (0.162466, (0.065771, 0.223679, 1.0, 1.0)),
                    (0.206899, (0.045999, 0.35081, 1.0, 1.0)),
                    (0.244976, (0.022557, 0.466392, 0.88704, 1.0)),
                    (0.299272, (0.005519, 0.679385, 0.591104, 1.0)),
                    (0.365893, (0.011259, 0.873389, 0.350493, 1.0)),
                    (0.426372, (0.091633, 0.99648, 0.129745, 1.0)),
                    (0.485072, (0.320069, 1.0, 0.045338, 1.0)),
                    (0.540405, (0.496488, 0.937007, 0.029813, 1.0)),
                    (0.60606, (0.791286, 0.689735, 0.040792, 1.0)),
                    (0.670056, (1.0, 0.489445, 0.045067, 1.0)),
                    (0.738035, (0.993354, 0.233808, 0.016193, 1.0)),
                    (0.809221, (0.879231, 0.070827, 0.003277, 1.0)),
                    (0.881656, (0.599744, 0.021882, 0.001288, 1.0)),
                    (0.93396, (0.453736, 0.007898, 0.0, 1.0)),
                    (1.0, (0.147729, 0.0, 0.000856, 1.0)),
                ),
                color_interpolation="B_SPLINE",
            )
        with g.Frame("nipy_spectral"):
            _string_81 = g.String(
                string="nipy_spectral: Miscellaneous. 32 stops fitted to matplotlib 3.11.2, linear interpolation, max error 5.0/255 in sRGB. Detail finer than the stops can hold is lost."
            )
            color_ramp_81 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.00011, 0.0, 0.0, 1.0)),
                    (0.010494, (0.007957, 0.000001, 0.009589, 1.0)),
                    (0.026408, (0.044638, 0.0, 0.057909, 1.0)),
                    (0.052157, (0.185571, 0.0, 0.247246, 1.0)),
                    (0.098309, (0.245303, 0.0, 0.317596, 1.0)),
                    (0.122412, (0.06469, 0.0, 0.350884, 1.0)),
                    (0.13639, (0.016702, 0.0, 0.384874, 1.0)),
                    (0.147915, (0.000406, 0.0, 0.387067, 1.0)),
                    (0.201584, (0.0, 0.000027, 0.722928, 1.0)),
                    (0.21417, (0.0, 0.013717, 0.723021, 1.0)),
                    (0.229422, (0.0, 0.056929, 0.722305, 1.0)),
                    (0.252267, (0.0, 0.189712, 0.724328, 1.0)),
                    (0.299003, (0.0, 0.312282, 0.721346, 1.0)),
                    (0.347866, (0.0, 0.401126, 0.406807, 1.0)),
                    (0.397987, (0.0, 0.405854, 0.249171, 1.0)),
                    (0.424165, (0.0, 0.354798, 0.054978, 1.0)),
                    (0.439136, (0.0, 0.346099, 0.011119, 1.0)),
                    (0.450181, (0.0, 0.297923, 0.0, 1.0)),
                    (0.601954, (0.000358, 0.968608, 0.0, 1.0)),
                    (0.610695, (0.016912, 1.0, 0.0, 1.0)),
                    (0.626428, (0.11015, 0.998028, 0.0, 1.0)),
                    (0.654642, (0.538005, 1.0, 0.0, 1.0)),
                    (0.716987, (0.938205, 0.772505, 0.0, 1.0)),
                    (0.795894, (1.0, 0.347043, 0.0, 1.0)),
                    (0.821688, (1.0, 0.084563, 0.0, 1.0)),
                    (0.838894, (1.0, 0.013816, 0.0, 1.0)),
                    (0.849654, (0.955246, 0.0, 0.0, 1.0)),
                    (0.950874, (0.567318, 0.0, 0.0, 1.0)),
                    (0.961371, (0.627362, 0.023662, 0.023662, 1.0)),
                    (0.973217, (0.590887, 0.107973, 0.107973, 1.0)),
                    (0.986449, (0.611965, 0.287147, 0.287147, 1.0)),
                    (1.0, (0.599039, 0.597632, 0.597632, 1.0)),
                ),
            )
        with g.Frame("gist_ncar"):
            _string_82 = g.String(
                string="gist_ncar: Miscellaneous. 32 stops fitted to matplotlib 3.11.2, linear interpolation, max error 6.4/255 in sRGB. Detail finer than the stops can hold is lost."
            )
            color_ramp_82 = g.ColorRamp(
                fac=switch,
                items=(
                    (0.0, (0.0, 0.0, 0.214185, 1.0)),
                    (0.015276, (0.0, 0.009725, 0.099006, 1.0)),
                    (0.037153, (0.0, 0.058455, 0.015291, 1.0)),
                    (0.053597, (0.0, 0.112783, 0.000291, 1.0)),
                    (0.069094, (0.0, 0.046181, 0.070692, 1.0)),
                    (0.087319, (0.0, 0.011005, 0.324181, 1.0)),
                    (0.108109, (0.0, 0.0, 0.943162, 1.0)),
                    (0.115334, (0.0, 0.013508, 1.0, 1.0)),
                    (0.132572, (0.0, 0.102528, 1.0, 1.0)),
                    (0.163174, (0.0, 0.545048, 0.999991, 1.0)),
                    (0.205737, (0.0, 0.902867, 0.995478, 1.0)),
                    (0.217614, (0.0, 1.0, 0.780889, 1.0)),
                    (0.284434, (0.0, 0.952135, 0.113754, 1.0)),
                    (0.305209, (0.0, 0.99986, 0.021292, 1.0)),
                    (0.319966, (0.003492, 0.989069, 0.0, 1.0)),
                    (0.341487, (0.030439, 0.797162, 0.0, 1.0)),
                    (0.370666, (0.121379, 0.615846, 0.0, 1.0)),
                    (0.425908, (0.225244, 0.991455, 0.001345, 1.0)),
                    (0.474354, (0.472434, 1.0, 0.037377, 1.0)),
                    (0.526953, (0.962858, 0.972052, 0.0, 1.0)),
                    (0.630735, (1.0, 0.480213, 0.002824, 1.0)),
                    (0.662786, (0.999999, 0.182143, 0.00289, 1.0)),
                    (0.68888, (1.0, 0.047193, 0.0, 1.0)),
                    (0.738606, (1.0, 0.0, 0.0, 1.0)),
                    (0.749624, (1.0, 0.001031, 0.032989, 1.0)),
                    (0.759337, (1.0, 0.0, 0.118336, 1.0)),
                    (0.772432, (1.0, 0.0, 0.331108, 1.0)),
                    (0.794579, (0.980758, 0.000147, 0.984563, 1.0)),
                    (0.815916, (0.645959, 0.007647, 0.998704, 1.0)),
                    (0.849218, (0.338906, 0.029824, 0.988804, 1.0)),
                    (0.901951, (0.857466, 0.209475, 0.853262, 1.0)),
                    (1.0, (0.977879, 0.884764, 0.996627, 1.0)),
                ),
            )
        menu_switch = g.MenuSwitch.color(
            uniform,
            {
                "viridis": (
                    color_ramp.o.color,
                    "Perceptually uniform sequential colormap from matplotlib",
                ),
                "plasma": (
                    color_ramp_1.o.color,
                    "Perceptually uniform sequential colormap from matplotlib",
                ),
                "inferno": (
                    color_ramp_2.o.color,
                    "Perceptually uniform sequential colormap from matplotlib",
                ),
                "magma": (
                    color_ramp_3.o.color,
                    "Perceptually uniform sequential colormap from matplotlib",
                ),
                "cividis": (
                    color_ramp_4.o.color,
                    "Perceptually uniform sequential colormap from matplotlib",
                ),
            },
        )
        menu_switch_1 = g.MenuSwitch.color(
            sequential,
            {
                "Greys": (color_ramp_5.o.color, "Sequential colormap from matplotlib"),
                "Purples": (
                    color_ramp_6.o.color,
                    "Sequential colormap from matplotlib",
                ),
                "Blues": (color_ramp_7.o.color, "Sequential colormap from matplotlib"),
                "Greens": (color_ramp_8.o.color, "Sequential colormap from matplotlib"),
                "Oranges": (
                    color_ramp_9.o.color,
                    "Sequential colormap from matplotlib",
                ),
                "Reds": (color_ramp_10.o.color, "Sequential colormap from matplotlib"),
                "YlOrBr": (
                    color_ramp_11.o.color,
                    "Sequential colormap from matplotlib",
                ),
                "YlOrRd": (
                    color_ramp_12.o.color,
                    "Sequential colormap from matplotlib",
                ),
                "OrRd": (color_ramp_13.o.color, "Sequential colormap from matplotlib"),
                "PuRd": (color_ramp_14.o.color, "Sequential colormap from matplotlib"),
                "RdPu": (color_ramp_15.o.color, "Sequential colormap from matplotlib"),
                "BuPu": (color_ramp_16.o.color, "Sequential colormap from matplotlib"),
                "GnBu": (color_ramp_17.o.color, "Sequential colormap from matplotlib"),
                "PuBu": (color_ramp_18.o.color, "Sequential colormap from matplotlib"),
                "YlGnBu": (
                    color_ramp_19.o.color,
                    "Sequential colormap from matplotlib",
                ),
                "PuBuGn": (
                    color_ramp_20.o.color,
                    "Sequential colormap from matplotlib",
                ),
                "BuGn": (color_ramp_21.o.color, "Sequential colormap from matplotlib"),
                "YlGn": (color_ramp_22.o.color, "Sequential colormap from matplotlib"),
            },
        )
        menu_switch_2 = g.MenuSwitch.color(
            sequential_2,
            {
                "binary": (
                    color_ramp_23.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "gray": (
                    color_ramp_24.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "bone": (
                    color_ramp_25.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "pink": (
                    color_ramp_26.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "spring": (
                    color_ramp_27.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "summer": (
                    color_ramp_28.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "autumn": (
                    color_ramp_29.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "winter": (
                    color_ramp_30.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "cool": (
                    color_ramp_31.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "Wistia": (
                    color_ramp_32.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "hot": (
                    color_ramp_33.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "afmhot": (
                    color_ramp_34.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "gist_heat": (
                    color_ramp_35.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
                "copper": (
                    color_ramp_36.o.color,
                    "Sequential (2) colormap from matplotlib",
                ),
            },
        )
        menu_switch_3 = g.MenuSwitch.color(
            diverging,
            {
                "PiYG": (color_ramp_37.o.color, "Diverging colormap from matplotlib"),
                "PRGn": (color_ramp_38.o.color, "Diverging colormap from matplotlib"),
                "BrBG": (color_ramp_39.o.color, "Diverging colormap from matplotlib"),
                "PuOr": (color_ramp_40.o.color, "Diverging colormap from matplotlib"),
                "RdGy": (color_ramp_41.o.color, "Diverging colormap from matplotlib"),
                "RdBu": (color_ramp_42.o.color, "Diverging colormap from matplotlib"),
                "RdYlBu": (color_ramp_43.o.color, "Diverging colormap from matplotlib"),
                "RdYlGn": (color_ramp_44.o.color, "Diverging colormap from matplotlib"),
                "Spectral": (
                    color_ramp_45.o.color,
                    "Diverging colormap from matplotlib",
                ),
                "coolwarm": (
                    color_ramp_46.o.color,
                    "Diverging colormap from matplotlib",
                ),
                "bwr": (color_ramp_47.o.color, "Diverging colormap from matplotlib"),
                "seismic": (
                    color_ramp_48.o.color,
                    "Diverging colormap from matplotlib",
                ),
                "berlin": (color_ramp_49.o.color, "Diverging colormap from matplotlib"),
                "managua": (
                    color_ramp_50.o.color,
                    "Diverging colormap from matplotlib",
                ),
                "vanimo": (color_ramp_51.o.color, "Diverging colormap from matplotlib"),
            },
        )
        menu_switch_4 = g.MenuSwitch.color(
            cyclic,
            {
                "twilight": (color_ramp_52.o.color, "Cyclic colormap from matplotlib"),
                "twilight_shifted": (
                    color_ramp_53.o.color,
                    "Cyclic colormap from matplotlib",
                ),
                "hsv": (color_ramp_54.o.color, "Cyclic colormap from matplotlib"),
            },
        )
        menu_switch_5 = g.MenuSwitch.color(
            qualitative,
            {
                "Pastel1": (
                    color_ramp_55.o.color,
                    "Qualitative colormap from matplotlib",
                ),
                "Pastel2": (
                    color_ramp_56.o.color,
                    "Qualitative colormap from matplotlib",
                ),
                "Paired": (
                    color_ramp_57.o.color,
                    "Qualitative colormap from matplotlib",
                ),
                "Accent": (
                    color_ramp_58.o.color,
                    "Qualitative colormap from matplotlib",
                ),
                "Dark2": (
                    color_ramp_59.o.color,
                    "Qualitative colormap from matplotlib",
                ),
                "Set1": (color_ramp_60.o.color, "Qualitative colormap from matplotlib"),
                "Set2": (color_ramp_61.o.color, "Qualitative colormap from matplotlib"),
                "Set3": (color_ramp_62.o.color, "Qualitative colormap from matplotlib"),
                "tab10": (
                    color_ramp_63.o.color,
                    "Qualitative colormap from matplotlib",
                ),
                "tab20": (
                    color_ramp_64.o.color,
                    "Qualitative colormap from matplotlib",
                ),
                "tab20b": (
                    color_ramp_65.o.color,
                    "Qualitative colormap from matplotlib",
                ),
                "tab20c": (
                    color_ramp_66.o.color,
                    "Qualitative colormap from matplotlib",
                ),
                "okabe_ito": (
                    color_ramp_67.o.color,
                    "Qualitative colormap from matplotlib",
                ),
            },
        )
        menu_switch_6 = g.MenuSwitch.color(
            miscellaneous,
            {
                "ocean": (
                    color_ramp_68.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "gist_earth": (
                    color_ramp_69.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "terrain": (
                    color_ramp_70.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "gist_stern": (
                    color_ramp_71.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "gnuplot": (
                    color_ramp_72.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "gnuplot2": (
                    color_ramp_73.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "CMRmap": (
                    color_ramp_74.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "cubehelix": (
                    color_ramp_75.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "brg": (
                    color_ramp_76.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "gist_rainbow": (
                    color_ramp_77.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "rainbow": (
                    color_ramp_78.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "jet": (
                    color_ramp_79.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "turbo": (
                    color_ramp_80.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "nipy_spectral": (
                    color_ramp_81.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
                "gist_ncar": (
                    color_ramp_82.o.color,
                    "Miscellaneous colormap from matplotlib",
                ),
            },
        )
        (
            g.MenuSwitch.color(
                category,
                {
                    "Uniform": (
                        menu_switch.o.output,
                        "Perceptually uniform sequential",
                    ),
                    "Sequential": (menu_switch_1.o.output, "Sequential"),
                    "Sequential 2": (menu_switch_2.o.output, "Sequential (2)"),
                    "Diverging": (menu_switch_3.o.output, "Diverging"),
                    "Cyclic": (menu_switch_4.o.output, "Cyclic"),
                    "Qualitative": (menu_switch_5.o.output, "Qualitative"),
                    "Miscellaneous": (menu_switch_6.o.output, "Miscellaneous"),
                },
            )
            >> color
        )

        category.default_value = "Uniform"
        uniform.default_value = "viridis"
        sequential.default_value = "Greys"
        sequential_2.default_value = "binary"
        diverging.default_value = "PiYG"
        cyclic.default_value = "twilight"
        qualitative.default_value = "Pastel1"
        miscellaneous.default_value = "ocean"


ASSET = ColorMatplotlib

ASSET_METADATA = {
    "description": "Colour a value from 0 to 1 with any of matplotlib's colormaps",
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
