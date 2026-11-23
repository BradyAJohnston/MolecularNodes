# Node-group asset "Color Sample Grid" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    FloatSocket,
    ObjectSocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputColor, InputFloat, InputObject, InputString


class ColorSampleGrid(AssetGeometryGroup):
    """
    Color Sample Grid

    Parameters
    ----------
    object : InputObject
        Object holding the grid to sample: a density loaded with Molecular Nodes, styled or not, or any volume object
    grid_name : InputString
        Name of the grid to sample: the Name of the density's Import Volume node (`density` by default), or a grid of the object's own volume
    min : InputFloat
        Grid value mapped to color A, in the grid's own units (kT/e for an APBS potential)
    max : InputFloat
        Grid value mapped to color B, in the grid's own units
    a : InputColor
        Color for values at or below Min
    intermediate : InputColor
        Color for the value halfway between Min and Max. Colors blend linearly from A to Intermediate to B
    b : InputColor
        Color for values at or above Max

    Inputs
    ------
    i.object : ObjectSocket
        Object holding the grid to sample: a density loaded with Molecular Nodes, styled or not, or any volume object
    i.grid_name : StringSocket
        Name of the grid to sample: the Name of the density's Import Volume node (`density` by default), or a grid of the object's own volume
    i.min : FloatSocket
        Grid value mapped to color A, in the grid's own units (kT/e for an APBS potential)
    i.max : FloatSocket
        Grid value mapped to color B, in the grid's own units
    i.a : ColorSocket
        Color for values at or below Min
    i.intermediate : ColorSocket
        Color for the value halfway between Min and Max. Colors blend linearly from A to Intermediate to B
    i.b : ColorSocket
        Color for values at or above Max

    Outputs
    -------
    o.color : ColorSocket
        Color mapped from the grid value at each point
    o.value : FloatSocket
        Grid value at each point, trilinearly interpolated. Zero outside the grid
    """

    _name = "Color Sample Grid"
    _asset_name = "Color Sample Grid"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        object: ObjectSocket
        """Object holding the grid to sample: a density loaded with Molecular Nodes, styled or not, or any volume object"""
        grid_name: StringSocket
        """Name of the grid to sample: the Name of the density's Import Volume node (`density` by default), or a grid of the object's own volume"""
        min: FloatSocket
        """Grid value mapped to color A, in the grid's own units (kT/e for an APBS potential)"""
        max: FloatSocket
        """Grid value mapped to color B, in the grid's own units"""
        a: ColorSocket
        """Color for values at or below Min"""
        intermediate: ColorSocket
        """Color for the value halfway between Min and Max. Colors blend linearly from A to Intermediate to B"""
        b: ColorSocket
        """Color for values at or above Max"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """Color mapped from the grid value at each point"""
        value: FloatSocket
        """Grid value at each point, trilinearly interpolated. Zero outside the grid"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        object: InputObject = None,
        grid_name: InputString = "density",
        min: InputFloat = -5.0,
        max: InputFloat = 5.0,
        a: InputColor = None,
        intermediate: InputColor = None,
        b: InputColor = None,
    ):
        super().__init__(
            **{
                "Object": object,
                "Grid Name": grid_name,
                "Min": min,
                "Max": max,
                "A": a,
                "Intermediate": intermediate,
                "B": b,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        object = tree.inputs.object(
            "Object",
            description="Object holding the grid to sample: a density loaded with Molecular Nodes, styled or not, or any volume object",
        )
        grid_name = tree.inputs.string(
            "Grid Name",
            "density",
            description="Name of the grid to sample: the Name of the density's Import Volume node (`density` by default), or a grid of the object's own volume",
        )
        min = tree.inputs.float(
            "Min",
            -5.0,
            description="Grid value mapped to color A, in the grid's own units (kT/e for an APBS potential)",
        )
        max = tree.inputs.float(
            "Max",
            5.0,
            description="Grid value mapped to color B, in the grid's own units",
        )
        with tree.inputs.panel("Color"):
            a = tree.inputs.color(
                "A",
                (0.6, 0.02, 0.02, 1.0),
                description="Color for values at or below Min",
            )
            intermediate = tree.inputs.color(
                "Intermediate",
                (0.9, 0.9, 0.9, 1.0),
                description="Color for the value halfway between Min and Max. Colors blend linearly from A to Intermediate to B",
            )
            b = tree.inputs.color(
                "B",
                (0.02, 0.08, 0.6, 1.0),
                description="Color for values at or above Max",
            )
        color = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 0.0),
            description="Color mapped from the grid value at each point",
        )
        value = tree.outputs.float(
            "Value",
            description="Grid value at each point, trilinearly interpolated. Zero outside the grid",
        )

        object_info = g.ObjectInfo(object=object, transform_space="RELATIVE")
        with g.Frame("A density keeps its grids in its bundle as 'MN Grid <name>'"):
            get_bundle_item = g.GetBundleItem.float(
                object_info.o.geometry.bundle(),
                g.FormatString("MN Grid {name}", items={"name": grid_name}),
                structure_type="GRID",
            )
            _string = g.String(
                "Import Volume stores its grid in the bundle of its output as 'MN Grid <name>', and density styles pass it on to the mesh they make. Sample that grid if it is there, otherwise the grid of that name in the object's own volume."
            )
            switch = get_bundle_item.o.exists.switch.float(
                g.GetNamedGrid.float(object_info.o.geometry, grid_name).o.grid,
                get_bundle_item.o.item,
            )
        sample_grid = g.SampleGrid.float(switch, g.Position())
        map_range = sample_grid.o.value.map_range(min, max)
        (
            (map_range * 2.0 - 1.0).mix.color(
                (map_range * 2.0).mix.color(a, intermediate), b
            )
            >> color
        )

        sample_grid >> value


ASSET = ColorSampleGrid

ASSET_METADATA = {
    "description": "Color points by sampling a float grid from a volume object, such as an electrostatic potential on a molecular surface",
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
