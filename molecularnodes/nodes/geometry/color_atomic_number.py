# Node-group asset 'Color Atomic Number' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputInteger
from .atomic_number import AtomicNumber
from .color import Color


class ColorAtomicNumber(AssetGeometryGroup):
    """
    Color Atomic Number

    Parameters
    ----------
    atomic_number : InputInteger
        The `atomic_number` of to use the selected color for
    color : InputColor
        The color to use for the specified `atomic_number`

    Inputs
    ------
    i.atomic_number : IntegerSocket
        The `atomic_number` of to use the selected color for
    i.color : ColorSocket
        The color to use for the specified `atomic_number`

    Outputs
    -------
    o.color : ColorSocket
        The generated color based on the node inputs
    """

    _name = "Color Atomic Number"
    _asset_name = "Color Atomic Number"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_atomic_number"}

    class _Inputs(SocketAccessor):
        atomic_number: IntegerSocket
        """The `atomic_number` of to use the selected color for"""
        color: ColorSocket
        """The color to use for the specified `atomic_number`"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The generated color based on the node inputs"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atomic_number: InputInteger = 6,
        color: InputColor = None,
    ):
        super().__init__(**{"atomic_number": atomic_number, "Color": color})

    def _build_group(self, tree):
        atomic_number = tree.inputs.integer(
            "atomic_number",
            6,
            description="The `atomic_number` of to use the selected color for",
            min_value=1,
            max_value=140,
        )
        color = tree.inputs.color(
            "Color",
            (0.8, 0.8, 0.8, 1.0),
            description="The color to use for the specified `atomic_number`",
        )
        color_1 = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 0.0),
            description="The generated color based on the node inputs",
        )

        (
            g.Compare.integer.equal(
                AtomicNumber(), atomic_number
            ).o.result.switch.color(Color(), color)
            >> color_1
        )


ASSET = ColorAtomicNumber

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
