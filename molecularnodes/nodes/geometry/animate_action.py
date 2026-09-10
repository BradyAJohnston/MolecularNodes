# Node-group asset 'Animate Action' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .between_integer import BetweenInteger


class AnimateAction(AssetGeometryGroup):
    """
    Animate Action

    Parameters
    ----------
    start : InputInteger
        Start
    length : InputInteger
        Length

    Inputs
    ------
    i.start : IntegerSocket
        Start
    i.length : IntegerSocket
        Length

    Outputs
    -------
    o.factor : FloatSocket
        Factor
    o.active : BooleanSocket
        Whether the input `Value` is between (and including) the lower and upper bounds
    o.stop : IntegerSocket
        Stop
    """

    _name = "Animate Action"
    _asset_name = "Animate Action"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        start: IntegerSocket
        """Start"""
        length: IntegerSocket
        """Length"""

    class _Outputs(SocketAccessor):
        factor: FloatSocket
        """Factor"""
        active: BooleanSocket
        """Whether the input `Value` is between (and including) the lower and upper bounds"""
        stop: IntegerSocket
        """Stop"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        start: InputInteger = 1,
        length: InputInteger = 100,
    ):
        super().__init__(**{"Start": start, "Length": length})

    def _build_group(self, tree):
        start = tree.inputs.integer("Start", 1)
        length = tree.inputs.integer("Length", 100, min_value=0)
        factor = tree.outputs.float("Factor")
        active = tree.outputs.boolean(
            "Active",
            description="Whether the input `Value` is between (and including) the lower and upper bounds",
        )
        stop = tree.outputs.integer("Stop")

        integer_math = start + length
        scene_time = g.SceneTime()
        (
            BetweenInteger(value=scene_time.o.frame, lower=start, upper=integer_math)
            >> active
        )
        scene_time.o.frame.map_range(start, integer_math) >> factor

        integer_math >> stop


ASSET = AnimateAction

ASSET_METADATA = {
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
