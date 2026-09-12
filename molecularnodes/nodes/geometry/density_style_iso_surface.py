# Node-group asset "Density Style ISO Surface" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputColor,
    InputFloat,
    InputGeometry,
    InputMaterial,
    InputVector,
)
from .between_vector import BetweenVector


class DensityStyleISOSurface(AssetGeometryGroup):
    """
    Density Style ISO Surface

    Parameters
    ----------
    volume : InputGeometry
        Input geometry
    visible : InputBoolean
        Visibility of style
    threshold : InputFloat
        ISO value
    show_contours : InputBoolean
        Whether to show surface contours
    only_contours : InputBoolean
        Only show contour edges
    contour_thickness : InputFloat
        Thickness of the contour edges
    contour_color : InputColor
        Color of contour edges
    slice_width : InputVector
        1
    slice_center : InputVector
        Slice Center
    negative_color : InputColor
        Color for negative ISO values
    positive_color : InputColor
        Color for positive ISO values
    shade_smooth : InputBoolean
        Use smooth shading for surface
    material : InputMaterial
        Material to use for this surface

    Inputs
    ------
    i.volume : GeometrySocket
        Input geometry
    i.visible : BooleanSocket
        Visibility of style
    i.threshold : FloatSocket
        ISO value
    i.show_contours : BooleanSocket
        Whether to show surface contours
    i.only_contours : BooleanSocket
        Only show contour edges
    i.contour_thickness : FloatSocket
        Thickness of the contour edges
    i.contour_color : ColorSocket
        Color of contour edges
    i.slice_width : VectorSocket
        1
    i.slice_center : VectorSocket
        Slice Center
    i.negative_color : ColorSocket
        Color for negative ISO values
    i.positive_color : ColorSocket
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
        visible: BooleanSocket
        """Visibility of style"""
        threshold: FloatSocket
        """ISO value"""
        show_contours: BooleanSocket
        """Whether to show surface contours"""
        only_contours: BooleanSocket
        """Only show contour edges"""
        contour_thickness: FloatSocket
        """Thickness of the contour edges"""
        contour_color: ColorSocket
        """Color of contour edges"""
        slice_width: VectorSocket
        """1"""
        slice_center: VectorSocket
        """Slice Center"""
        negative_color: ColorSocket
        """Color for negative ISO values"""
        positive_color: ColorSocket
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
        visible: InputBoolean = True,
        threshold: InputFloat = 0.00959343,
        show_contours: InputBoolean = False,
        only_contours: InputBoolean = False,
        contour_thickness: InputFloat = 0.1,
        contour_color: InputColor = None,
        slice_width: InputVector = None,
        slice_center: InputVector = None,
        negative_color: InputColor = None,
        positive_color: InputColor = None,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Volume": volume,
                "Visible": visible,
                "Threshold": threshold,
                "Show Contours": show_contours,
                "Only Contours": only_contours,
                "Contour Thickness": contour_thickness,
                "Contour Color": contour_color,
                "Slice Width": slice_width,
                "Slice Center": slice_center,
                "Negative Color": negative_color,
                "Positive Color": positive_color,
                "Shade Smooth": shade_smooth,
                "Material": material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        volume = tree.inputs.geometry("Volume", description="Input geometry")
        visible = tree.inputs.boolean(
            "Visible", True, description="Visibility of style"
        )
        threshold = tree.inputs.float(
            "Threshold",
            0.00959343,
            description="ISO value",
            min_value=0.0,
            max_value=1.0,
        )
        with tree.inputs.panel("Contours", default_closed=True):
            show_contours = tree.inputs.boolean(
                "Show Contours", False, description="Whether to show surface contours"
            )
            only_contours = tree.inputs.boolean(
                "Only Contours", False, description="Only show contour edges"
            )
            contour_thickness = tree.inputs.float(
                "Contour Thickness",
                0.1,
                description="Thickness of the contour edges",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            contour_color = tree.inputs.color(
                "Contour Color",
                (0.0, 0.0, 0.0, 1.0),
                description="Color of contour edges",
            )
        with tree.inputs.panel("Slice"):
            slice_width = tree.inputs.vector(
                "Slice Width",
                (0.5, 0.5, 0.5),
                description="1",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            slice_center = tree.inputs.vector(
                "Slice Center",
                (0.5, 0.5, 0.5),
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
        with tree.inputs.panel("Material", default_closed=True):
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
            material = tree.inputs.material(
                "Material", description="Material to use for this surface"
            )
        geometry = tree.outputs.geometry(
            "Geometry", description="ISO surface geometry output"
        )

        bounding_box = g.BoundingBox(geometry=volume)
        math_1 = contour_thickness * 0.001
        switch = visible.switch.geometry(true=volume)
        group = BetweenVector(
            value=g.Position().o.position.map_range(
                bounding_box.o.min, bounding_box.o.max
            ),
            lower=slice_center - slice_width,
            upper=slice_center + slice_width,
        )
        set_material = (
            g.VolumeToMesh(volume=switch, threshold=threshold * -1.0, voxel_size=0.3)
            >> g.StoreNamedAttribute.point.color(name="Color", value=negative_color)
            >> g.SetMaterial(material=material)
        )
        set_material_1 = (
            g.VolumeToMesh(volume=switch, threshold=threshold, voxel_size=0.3)
            >> g.StoreNamedAttribute.point.color(name="Color", value=positive_color)
            >> g.SetMaterial(material=material)
        )
        separate_geometry = (
            g.JoinGeometry(geometry=(set_material_1, set_material))
            >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
            >> g.SeparateGeometry.point(selection=group)
        )
        set_material_2 = (
            g.MeshToCurve(mesh=separate_geometry.o.selection, selection=show_contours)
            >> g.CurveToMesh(profile_curve=g.Quadrilateral(width=math_1, height=math_1))
            >> g.StoreNamedAttribute.edge.color(name="Color", value=contour_color)
            >> g.SetMaterial(material=material)
        )
        join_geometry = g.JoinGeometry(
            geometry=(
                only_contours.switch.geometry(separate_geometry.o.selection),
                set_material_2,
            )
        )

        join_geometry >> geometry


ASSET = DensityStyleISOSurface

ASSET_METADATA = {
    "catalog_id": "35b7cc56-45fd-4113-8e88-3c7387edb2d8",
}
