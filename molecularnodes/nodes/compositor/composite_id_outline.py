# Node-group asset "Composite ID Outline" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetCompositorGroup,
    CustomCompositorGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat


class IDDifferenceSample(CustomCompositorGroup):
    _name = "ID Difference Sample"
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Whether the ID at a pixel offset differs from the centre pixel's ID by more than Min Difference"
    }

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


class CompositeIDOutline(AssetCompositorGroup):
    """
    Outline opacity between regions of an ID pass after Goodsell's Illustrate: counts, over the 5x5 neighbourhood, the pixels whose ID differs from the centre by more than Min Difference and ramps the count from Low to High into an opacity. Feed it a float AOV such as chain_id or res_id (canvas.add_aov and mn.material.add_aov) for chain or residue outlines

    Parameters
    ----------
    id : InputFloat
        A float pass identifying regions, such as a chain_id or res_id AOV from the Render Layers node
    min_difference : InputFloat
        Difference in ID above which two pixels belong to different regions. 0.5 separates integer IDs; raise it to outline every few residues of a res_id pass instead of each one
    low : InputFloat
        Number of differing neighbours, of 24, at which the line starts to appear
    high : InputFloat
        Number of differing neighbours, of 24, at which the line is fully opaque

    Inputs
    ------
    i.id : FloatSocket
        A float pass identifying regions, such as a chain_id or res_id AOV from the Render Layers node
    i.min_difference : FloatSocket
        Difference in ID above which two pixels belong to different regions. 0.5 separates integer IDs; raise it to outline every few residues of a res_id pass instead of each one
    i.low : FloatSocket
        Number of differing neighbours, of 24, at which the line starts to appear
    i.high : FloatSocket
        Number of differing neighbours, of 24, at which the line is fully opaque

    Outputs
    -------
    o.opacity : FloatSocket
        Opacity
    """

    _name = "Composite ID Outline"
    _asset_name = "Composite ID Outline"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Outline opacity between regions of an ID pass after Goodsell's Illustrate: counts, over the 5x5 neighbourhood, the pixels whose ID differs from the centre by more than Min Difference and ramps the count from Low to High into an opacity. Feed it a float AOV such as chain_id or res_id (canvas.add_aov and mn.material.add_aov) for chain or residue outlines"
    }

    class _Inputs(SocketAccessor):
        id: FloatSocket
        """A float pass identifying regions, such as a chain_id or res_id AOV from the Render Layers node"""
        min_difference: FloatSocket
        """Difference in ID above which two pixels belong to different regions. 0.5 separates integer IDs; raise it to outline every few residues of a res_id pass instead of each one"""
        low: FloatSocket
        """Number of differing neighbours, of 24, at which the line starts to appear"""
        high: FloatSocket
        """Number of differing neighbours, of 24, at which the line is fully opaque"""

    class _Outputs(SocketAccessor):
        opacity: FloatSocket
        """Opacity"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        id: InputFloat = 0.0,
        min_difference: InputFloat = 0.5,
        low: InputFloat = 3.0,
        high: InputFloat = 10.0,
    ):
        super().__init__(
            **{"ID": id, "Min Difference": min_difference, "Low": low, "High": high}
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        id = tree.inputs.float(
            "ID",
            0.0,
            description="A float pass identifying regions, such as a chain_id or res_id AOV from the Render Layers node",
            hide_value=True,
        )
        min_difference = tree.inputs.float(
            "Min Difference",
            0.5,
            description="Difference in ID above which two pixels belong to different regions. 0.5 separates integer IDs; raise it to outline every few residues of a res_id pass instead of each one",
            min_value=0.0,
            max_value=1_000_000.0,
        )
        low = tree.inputs.float(
            "Low",
            3.0,
            description="Number of differing neighbours, of 24, at which the line starts to appear",
            min_value=0.0,
            max_value=24.0,
        )
        high = tree.inputs.float(
            "High",
            10.0,
            description="Number of differing neighbours, of 24, at which the line is fully opaque",
            min_value=0.0,
            max_value=24.0,
        )
        opacity = tree.outputs.float("Opacity")

        with c.Frame("Count differing neighbours"):
            math_1 = IDDifferenceSample(
                ID=id, X=-2.0, Y=-2.0, **{"Min Difference": min_difference}
            ).o.differs + IDDifferenceSample(
                ID=id, X=-1.0, Y=-2.0, **{"Min Difference": min_difference}
            )
            math_2 = (
                math_1
                + IDDifferenceSample(
                    ID=id, Y=-2.0, **{"Min Difference": min_difference}
                )
                + IDDifferenceSample(
                    ID=id, X=1.0, Y=-2.0, **{"Min Difference": min_difference}
                )
            )
            math_3 = (
                math_2
                + IDDifferenceSample(
                    ID=id, X=2.0, Y=-2.0, **{"Min Difference": min_difference}
                )
                + IDDifferenceSample(
                    ID=id, X=-2.0, Y=-1.0, **{"Min Difference": min_difference}
                )
            )
            math_4 = (
                math_3
                + IDDifferenceSample(
                    ID=id, X=-1.0, Y=-1.0, **{"Min Difference": min_difference}
                )
                + IDDifferenceSample(
                    ID=id, Y=-1.0, **{"Min Difference": min_difference}
                )
            )
            math_5 = (
                math_4
                + IDDifferenceSample(
                    ID=id, X=1.0, Y=-1.0, **{"Min Difference": min_difference}
                )
                + IDDifferenceSample(
                    ID=id, X=2.0, Y=-1.0, **{"Min Difference": min_difference}
                )
            )
            math_6 = (
                math_5
                + IDDifferenceSample(
                    ID=id, X=-2.0, **{"Min Difference": min_difference}
                )
                + IDDifferenceSample(
                    ID=id, X=-1.0, **{"Min Difference": min_difference}
                )
            )
            math_7 = (
                math_6
                + IDDifferenceSample(ID=id, X=1.0, **{"Min Difference": min_difference})
                + IDDifferenceSample(ID=id, X=2.0, **{"Min Difference": min_difference})
            )
            math_8 = (
                math_7
                + IDDifferenceSample(
                    ID=id, X=-2.0, Y=1.0, **{"Min Difference": min_difference}
                )
                + IDDifferenceSample(
                    ID=id, X=-1.0, Y=1.0, **{"Min Difference": min_difference}
                )
            )
            math_9 = (
                math_8
                + IDDifferenceSample(ID=id, Y=1.0, **{"Min Difference": min_difference})
                + IDDifferenceSample(
                    ID=id, X=1.0, Y=1.0, **{"Min Difference": min_difference}
                )
            )
            math_10 = (
                math_9
                + IDDifferenceSample(
                    ID=id, X=2.0, Y=1.0, **{"Min Difference": min_difference}
                )
                + IDDifferenceSample(
                    ID=id, X=-2.0, Y=2.0, **{"Min Difference": min_difference}
                )
            )
            math_11 = (
                math_10
                + IDDifferenceSample(
                    ID=id, X=-1.0, Y=2.0, **{"Min Difference": min_difference}
                )
                + IDDifferenceSample(ID=id, Y=2.0, **{"Min Difference": min_difference})
            )
            math_12 = (
                math_11
                + IDDifferenceSample(
                    ID=id, X=1.0, Y=2.0, **{"Min Difference": min_difference}
                )
                + IDDifferenceSample(
                    ID=id, X=2.0, Y=2.0, **{"Min Difference": min_difference}
                )
            )
            _string = g.String(
                string="The 24 pixels of the 5x5 window around each pixel are compared with its ID; the number of pixels in another region is the line strength."
            )
        math_12.map_range(low, high) >> opacity


ASSET = CompositeIDOutline

ASSET_METADATA = {
    "description": "Outlines between regions of an ID pass, for chain and residue outlines",
    "catalog_id": "441e6ca5-e514-4e77-a3cd-25fc1a2e08ae",
}
