# Node group "ID Difference Sample" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy import geometry as g
from nodebpy.builder import CustomCompositorGroup, FloatSocket, SocketAccessor
from nodebpy.types import InputFloat


class IDDifferenceSample(CustomCompositorGroup):
    """
    Whether the ID at a pixel offset differs from the centre pixel's ID by more than Min Difference

    Parameters
    ----------
    id : InputFloat
        ID pass
    x : InputFloat
        Offset of the sample in pixels
    y : InputFloat
        Offset of the sample in pixels
    min_difference : InputFloat
        Difference in ID above which the sample counts as another region

    Inputs
    ------
    i.id : FloatSocket
        ID pass
    i.x : FloatSocket
        Offset of the sample in pixels
    i.y : FloatSocket
        Offset of the sample in pixels
    i.min_difference : FloatSocket
        Difference in ID above which the sample counts as another region

    Outputs
    -------
    o.differs : FloatSocket
        Differs
    """

    _name = "ID Difference Sample"
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Whether the ID at a pixel offset differs from the centre pixel's ID by more than Min Difference"
    }

    class _Inputs(SocketAccessor):
        id: FloatSocket
        """ID pass"""
        x: FloatSocket
        """Offset of the sample in pixels"""
        y: FloatSocket
        """Offset of the sample in pixels"""
        min_difference: FloatSocket
        """Difference in ID above which the sample counts as another region"""

    class _Outputs(SocketAccessor):
        differs: FloatSocket
        """Differs"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        id: InputFloat = 0.0,
        x: InputFloat = 0.0,
        y: InputFloat = 0.0,
        min_difference: InputFloat = 0.5,
    ):
        super().__init__(**{"ID": id, "X": x, "Y": y, "Min Difference": min_difference})

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        id = tree.inputs.float("ID", 0.0, description="ID pass", hide_value=True)
        x = tree.inputs.float("X", 0.0, description="Offset of the sample in pixels")
        y = tree.inputs.float("Y", 0.0, description="Offset of the sample in pixels")
        min_difference = tree.inputs.float(
            "Min Difference",
            0.5,
            description="Difference in ID above which the sample counts as another region",
        )
        differs = tree.outputs.float("Differs")

        translate = c.Translate(
            image=id,
            x=x,
            y=y,
            interpolation="Nearest",
            extension_x="Extend",
            extension_y="Extend",
        )
        math_1 = g.Math.greater_than(
            abs(g.Math.subtract(translate, id).o.value), min_difference
        )

        math_1 >> differs
