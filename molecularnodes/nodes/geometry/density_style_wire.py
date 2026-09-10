# Node-group asset 'Density Style Wire' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MaterialSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import (
    InputColor,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMaterial,
)


class DensityStyleWire(AssetGeometryGroup):
    """
    Density Style Wire

    Parameters
    ----------
    volume : InputGeometry
        Volume
    threshold : InputFloat
        Threshold
    hide_dust : InputFloat
        Hide Dust
    wire_radius : InputFloat
        Radius of the created wire (in relative nm)
    wire_resolution : InputInteger
        Wire Resolution
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
    i.hide_dust : FloatSocket
        Hide Dust
    i.wire_radius : FloatSocket
        Radius of the created wire (in relative nm)
    i.wire_resolution : IntegerSocket
        Wire Resolution
    i.color : ColorSocket
        Color
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        The generated geometry for the style node group
    """

    _name = "Density Style Wire"
    _asset_name = "Density Style Wire"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.style_density_wire"}

    class _Inputs(SocketAccessor):
        volume: GeometrySocket
        """Volume"""
        threshold: FloatSocket
        """Threshold"""
        hide_dust: FloatSocket
        """Hide Dust"""
        wire_radius: FloatSocket
        """Radius of the created wire (in relative nm)"""
        wire_resolution: IntegerSocket
        """Wire Resolution"""
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
        hide_dust: InputFloat = 20.0,
        wire_radius: InputFloat = 1.0,
        wire_resolution: InputInteger = 3,
        color: InputColor = None,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Volume": volume,
                "Threshold": threshold,
                "Hide Dust": hide_dust,
                "Wire Radius": wire_radius,
                "Wire Resolution": wire_resolution,
                "Color": color,
                "Material": material,
            }
        )

    def _build_group(self, tree):
        volume = tree.inputs.geometry("Volume")
        threshold = tree.inputs.float("Threshold", 0.8)
        hide_dust = tree.inputs.float(
            "Hide Dust", 20.0, min_value=-10_000.0, max_value=10_000.0
        )
        with tree.inputs.panel("Wire"):
            wire_radius = tree.inputs.float(
                "Wire Radius",
                1.0,
                description="Radius of the created wire (in relative nm)",
                min_value=0.0,
            )
            wire_resolution = tree.inputs.integer(
                "Wire Resolution", 3, min_value=3, max_value=32
            )
        with tree.inputs.panel("Material"):
            color = tree.inputs.color("Color", (0.10174983, 0.3931146, 0.10474136, 1.0))
            material = tree.inputs.material(
                "Material", description="Material to apply to the resulting geometry"
            )
        geometry = tree.outputs.geometry(
            "Geometry", description="The generated geometry for the style node group"
        )

        named_attribute = g.NamedAttribute.float("radius")
        (
            volume
            >> g.VolumeToMesh(threshold=threshold, voxel_size=0.3)
            >> g.DeleteGeometry.point(
                selection=g.FaceArea().o.area.point.total(g.MeshIsland().o.island_index)
                < hide_dust
            )
            >> g.StoreNamedAttribute.point.color(name="Color", value=color)
            >> g.MeshToCurve()
            >> g.CurveToMesh(
                profile_curve=g.CurveCircle(
                    resolution=wire_resolution, radius=wire_radius / 1000.0
                ),
                scale=named_attribute.o.exists.switch.float(
                    1.0, named_attribute.o.attribute
                ),
                fill_caps=True,
            )
            >> g.SetMaterial(material=material)
            >> geometry
        )


ASSET = DensityStyleWire

ASSET_METADATA = {
    "catalog_id": "35b7cc56-45fd-4113-8e88-3c7387edb2d8",
}
