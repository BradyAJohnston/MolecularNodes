# Node group '.Cleanup' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    GeometrySocket,
    MaterialSocket,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputGeometry, InputMaterial
from ..color import Color
from ..set_color import SetColor


class Cleanup(CustomGeometryGroup):
    """
    .Cleanup

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    color_source : InputGeometry
        Color Source
    material : InputMaterial
        Material to apply to the resulting geometry
    shade_smooth : InputBoolean
        Apply smooth shading to the created geometry

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.color_source : GeometrySocket
        Color Source
    i.material : MaterialSocket
        Material to apply to the resulting geometry
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = ".Cleanup"
    _tree_properties = {"node_tool_idname": "geometry._cleanup"}

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        color_source: GeometrySocket
        """Color Source"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""
        shade_smooth: BooleanSocket
        """Apply smooth shading to the created geometry"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        color_source: InputGeometry = None,
        material: InputMaterial = None,
        shade_smooth: InputBoolean = True,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Color Source": color_source,
                "Material": material,
                "Shade Smooth": shade_smooth,
            }
        )

    def _build_group(self, tree):
        geometry = tree.inputs.geometry("Geometry")
        color_source = tree.inputs.geometry("Color Source")
        material = tree.inputs.material(
            "Material", description="Material to apply to the resulting geometry"
        )
        shade_smooth = tree.inputs.boolean(
            "Shade Smooth",
            True,
            description="Apply smooth shading to the created geometry",
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        sample_index = g.SampleIndex(
            geometry=color_source,
            value=Color(),
            index=g.NamedAttribute.integer("tmp_idx").o.attribute,
            data_type="FLOAT_COLOR",
        )
        (
            SetColor(atoms=geometry, color=sample_index)
            >> g.RemoveNamedAttribute(pattern_mode="Wildcard", name="tmp_*")
            >> g.RemoveNamedAttribute(pattern_mode="Wildcard", name="backbone_*")
            >> g.SetMaterial(material=material)
            >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
            >> geometry_1
        )
