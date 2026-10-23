# Node-group asset "Density Style ISO Surface" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ClosureSocket,
    ColorSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    MaterialSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import (
    InputBoolean,
    InputClosure,
    InputColor,
    InputFloat,
    InputGeometry,
    InputMaterial,
)
from .set_color import SetColor


class IsoGridToMesh(CustomGeometryGroup):
    _name = "Iso Grid to Mesh"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        grid = tree.inputs.float("Grid", 0.0)
        process_grid = tree.inputs.closure("Process Grid")
        threshold = tree.inputs.float(
            "Threshold",
            0.00959343,
            description="ISO value. Surfaces are meshed at +Threshold and, for signed grids such as potentials or orbitals, at -Threshold",
        )
        material = tree.inputs.material(
            "Material", description="Material to use for this surface"
        )
        negative_color = tree.inputs.color(
            "Negative Color",
            (1.0, 0.0, 0.0, 1.0),
            description="Color for negative ISO values",
        )
        positive_color = tree.inputs.color(
            "Positive Color",
            (0.0, 0.0, 1.0, 1.0),
            description="Color for positive ISO values",
        )
        shade_smooth = tree.inputs.boolean(
            "Shade Smooth", True, description="Use smooth shading for surface"
        )
        geometry = tree.outputs.geometry("Geometry")

        closure_zone = g.ClosureZone()
        grid_1 = closure_zone.inputs.float("Grid")
        threshold_1 = closure_zone.inputs.float("Threshold")
        value = closure_zone.inputs.color("Value")
        negative = closure_zone.inputs.boolean("Negative")
        mesh = closure_zone.outputs.geometry("Mesh")
        evaluate_closure = g.EvaluateClosure(process_grid, define_signature=True)
        evaluate_closure.inputs.float(grid_1.output, "Grid")
        grid_2 = evaluate_closure.outputs.float("Grid")
        grid_to_mesh = g.GridToMesh(
            grid=grid_2.output,
            threshold=negative.output.switch.float(
                abs(threshold_1.output), abs(threshold_1.output) * -1.0
            ),
        )
        set_color = SetColor(
            atoms=g.FlipFaces(mesh=grid_to_mesh, selection=negative.output),
            color=value.output,
        )
        set_color >> mesh.input
        evaluate_closure_1 = g.EvaluateClosure(closure_zone.closure)
        evaluate_closure_1.inputs.float(grid, "Grid")
        evaluate_closure_1.inputs.float(threshold, "Threshold")
        evaluate_closure_1.inputs.color(negative_color, "Value")
        evaluate_closure_1.inputs.boolean(True, "Negative")
        mesh_1 = evaluate_closure_1.outputs.geometry("Mesh")
        evaluate_closure_2 = g.EvaluateClosure(closure_zone.closure)
        evaluate_closure_2.inputs.float(grid, "Grid")
        evaluate_closure_2.inputs.float(threshold, "Threshold")
        evaluate_closure_2.inputs.color(positive_color, "Value")
        evaluate_closure_2.inputs.boolean(name="Negative")
        mesh_2 = evaluate_closure_2.outputs.geometry("Mesh")
        (
            g.JoinGeometry((mesh_1.output, mesh_2.output))
            >> g.SetMaterial(material=material)
            >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
            >> g.JoinGeometry()
            >> geometry
        )


class GetGeometryGrid(CustomGeometryGroup):
    _name = "Get Geometry Grid"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to get the bundle of"
        )
        geometry_1 = tree.outputs.geometry(
            "Geometry", description="Geometry to get the bundle of"
        )
        bundle = tree.outputs.bundle("Bundle")
        grid = tree.outputs.float("Grid")

        get_bundle_item = g.GetBundleItem.float(
            geometry.bundle(), "MN/Density", structure_type="GRID"
        )
        (
            get_bundle_item.o.exists.switch.float(
                g.GetNamedGrid.float(geometry, "density").o.grid, get_bundle_item.o.item
            )
            >> grid
        )

        geometry >> geometry_1
        get_bundle_item >> bundle


class EvaluateOnGrid(CustomGeometryGroup):
    _name = "Evaluate on Grid"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to get the bundle of"
        )
        closure = tree.inputs.closure("Closure")
        geometry_1 = tree.outputs.geometry("Geometry")

        get_geometry_grid = GetGeometryGrid(Geometry=geometry)
        evaluate_closure = g.EvaluateClosure(closure)
        evaluate_closure.inputs.float(get_geometry_grid.o.grid, "Grid")
        geometry_2 = evaluate_closure.outputs.geometry("Geometry")
        (
            geometry_2.output.set_bundle(
                g.StoreBundleItem.float(
                    get_geometry_grid.o.bundle, "MN/Grid", get_geometry_grid.o.grid
                )
            )
            >> geometry_1
        )


class DensityStyleISOSurface(AssetGeometryGroup):
    """
    Density Style ISO Surface

    Parameters
    ----------
    volume : InputGeometry
        Input geometry
    threshold : InputFloat
        ISO value. Surfaces are meshed at +Threshold and, for signed grids such as potentials or orbitals, at -Threshold
    process_grid : InputClosure
        Extra operations on the Float grid before being turned into a mesh
    negative : InputColor
        Color for negative ISO values
    positive : InputColor
        Color for positive ISO values
    shade_smooth : InputBoolean
        Use smooth shading for surface
    material : InputMaterial
        Material to use for this surface

    Inputs
    ------
    i.volume : GeometrySocket
        Input geometry
    i.threshold : FloatSocket
        ISO value. Surfaces are meshed at +Threshold and, for signed grids such as potentials or orbitals, at -Threshold
    i.process_grid : ClosureSocket
        Extra operations on the Float grid before being turned into a mesh
    i.negative : ColorSocket
        Color for negative ISO values
    i.positive : ColorSocket
        Color for positive ISO values
    i.shade_smooth : BooleanSocket
        Use smooth shading for surface
    i.material : MaterialSocket
        Material to use for this surface

    Outputs
    -------
    o.geometry : GeometrySocket
        ISO surface geometry output
    """

    _name = "Density Style ISO Surface"
    _asset_name = "Density Style ISO Surface"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        volume: GeometrySocket
        """Input geometry"""
        threshold: FloatSocket
        """ISO value. Surfaces are meshed at +Threshold and, for signed grids such as potentials or orbitals, at -Threshold"""
        process_grid: ClosureSocket
        """Extra operations on the Float grid before being turned into a mesh"""
        negative: ColorSocket
        """Color for negative ISO values"""
        positive: ColorSocket
        """Color for positive ISO values"""
        shade_smooth: BooleanSocket
        """Use smooth shading for surface"""
        material: MaterialSocket
        """Material to use for this surface"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """ISO surface geometry output"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        volume: InputGeometry = None,
        threshold: InputFloat = 0.00959343,
        process_grid: InputClosure = None,
        negative: InputColor = None,
        positive: InputColor = None,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Volume": volume,
                "Threshold": threshold,
                "Process Grid": process_grid,
                "Negative": negative,
                "Positive": positive,
                "Shade Smooth": shade_smooth,
                "Material": material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        volume = tree.inputs.geometry("Volume", description="Input geometry")
        threshold = tree.inputs.float(
            "Threshold",
            0.00959343,
            description="ISO value. Surfaces are meshed at +Threshold and, for signed grids such as potentials or orbitals, at -Threshold",
        )
        with tree.inputs.panel("Surface", default_closed=True):
            process_grid = tree.inputs.closure(
                "Process Grid",
                description="Extra operations on the Float grid before being turned into a mesh",
            )
            negative = tree.inputs.color(
                "Negative",
                (1.0, 0.0, 0.0, 1.0),
                description="Color for negative ISO values",
            )
            positive = tree.inputs.color(
                "Positive",
                (0.0, 0.0, 1.0, 1.0),
                description="Color for positive ISO values",
            )
            shade_smooth = tree.inputs.boolean(
                "Shade Smooth", True, description="Use smooth shading for surface"
            )
            material = tree.inputs.material(
                "Material", description="Material to use for this surface"
            )
        geometry = tree.outputs.geometry(
            "Geometry", description="ISO surface geometry output"
        )

        closure_zone = g.ClosureZone()
        grid = closure_zone.inputs.float("Grid")
        geometry_1 = closure_zone.outputs.geometry("Geometry")
        iso_grid_to_mesh = IsoGridToMesh(
            Grid=grid.output,
            **{"Process Grid": process_grid},
            Threshold=threshold,
            Material=material,
            **{
                "Negative Color": negative,
                "Positive Color": positive,
                "Shade Smooth": shade_smooth,
            },
        )
        iso_grid_to_mesh >> geometry_1.input
        EvaluateOnGrid(Geometry=volume, Closure=closure_zone.closure) >> geometry


ASSET = DensityStyleISOSurface

ASSET_METADATA = {
    "catalog_id": "35b7cc56-45fd-4113-8e88-3c7387edb2d8",
}
