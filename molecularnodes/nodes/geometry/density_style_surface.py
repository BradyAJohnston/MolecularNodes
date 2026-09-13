# Node-group asset "Density Style Surface" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ColorSocket,
    FloatSocket,
    GeometrySocket,
    MaterialSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import (
    InputBoolean,
    InputColor,
    InputFloat,
    InputGeometry,
    InputMaterial,
)


class DensityStyleSurface(AssetGeometryGroup):
    """
    Density Style Surface

    Parameters
    ----------
    volume : InputGeometry
        Volume
    threshold : InputFloat
        Threshold
    shade_smooth : InputBoolean
        Apply smooth shading to the created geometry
    hide_dust : InputFloat
        Hide Dust
    color : InputColor
        Color
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.volume : GeometrySocket
        Volume
    i.threshold : FloatSocket
        Threshold
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry
    i.hide_dust : FloatSocket
        Hide Dust
    i.color : ColorSocket
        Color
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        The generated geometry for the style node group
    """

    _name = "Density Style Surface"
    _asset_name = "Density Style Surface"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.style_density_surface"}

    class _Inputs(SocketAccessor):
        volume: GeometrySocket
        """Volume"""
        threshold: FloatSocket
        """Threshold"""
        shade_smooth: BooleanSocket
        """Apply smooth shading to the created geometry"""
        hide_dust: FloatSocket
        """Hide Dust"""
        color: ColorSocket
        """Color"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """The generated geometry for the style node group"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        volume: InputGeometry = None,
        threshold: InputFloat = 0.8,
        shade_smooth: InputBoolean = True,
        hide_dust: InputFloat = 0.0,
        color: InputColor = None,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Volume": volume,
                "Threshold": threshold,
                "Shade Smooth": shade_smooth,
                "Hide Dust": hide_dust,
                "Color": color,
                "Material": material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        volume = tree.inputs.geometry("Volume")
        threshold = tree.inputs.float("Threshold", 0.8)
        shade_smooth = tree.inputs.boolean(
            "Shade Smooth",
            True,
            description="Apply smooth shading to the created geometry",
        )
        hide_dust = tree.inputs.float(
            "Hide Dust", 0.0, min_value=-10_000.0, max_value=10_000.0
        )
        with tree.inputs.panel("Material"):
            color = tree.inputs.color("Color", (0.199436, 0.509163, 0.13218, 1.0))
            material = tree.inputs.material(
                "Material", description="Material to apply to the resulting geometry"
            )
        geometry = tree.outputs.geometry(
            "Geometry", description="The generated geometry for the style node group"
        )

        (
            volume
            >> g.VolumeToMesh(threshold=threshold, voxel_size=0.3)
            >> g.DeleteGeometry.point(
                selection=g.FaceArea().o.area.point.total(g.MeshIsland().o.island_index)
                < hide_dust
            )
            >> g.StoreNamedAttribute.point.color(name="Color", value=color)
            >> g.SetMaterial(material=material)
            >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
            >> geometry
        )


ASSET = DensityStyleSurface

ASSET_METADATA = {
    "catalog_id": "35b7cc56-45fd-4113-8e88-3c7387edb2d8",
}
