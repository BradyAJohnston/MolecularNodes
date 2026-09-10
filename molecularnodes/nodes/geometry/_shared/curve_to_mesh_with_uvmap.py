# Node group 'Curve to Mesh with UVMap' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    GeometrySocket,
    IntegerSocket,
    MenuSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputGeometry, InputInteger, InputMenu
from ..fallback_geometry import FallbackGeometry
from .check_end_face_corner import CheckEndFaceCorner


class CurveToMeshWithUVMap(CustomGeometryGroup):
    """
    Curve to Mesh with UVMap

    Parameters
    ----------
    curve : InputGeometry
        Curve
    u_component : InputMenu | Literal["Factor", "Length"]
        U Component
    profile_resolution : InputInteger
        Profile Resolution
    fill_caps : InputBoolean
        Fill Caps

    Inputs
    ------
    i.curve : GeometrySocket
        Curve
    i.u_component : MenuSocket
        U Component
    i.profile_resolution : IntegerSocket
        Profile Resolution
    i.fill_caps : BooleanSocket
        Fill Caps

    Outputs
    -------
    o.mesh : GeometrySocket
        Mesh
    o.uv_map : VectorSocket
        uv_map
    """

    _name = "Curve to Mesh with UVMap"

    class _Inputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""
        u_component: MenuSocket
        """U Component"""
        profile_resolution: IntegerSocket
        """Profile Resolution"""
        fill_caps: BooleanSocket
        """Fill Caps"""

    class _Outputs(SocketAccessor):
        mesh: GeometrySocket
        """Mesh"""
        uv_map: VectorSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        curve: InputGeometry = None,
        u_component: InputMenu | Literal["Factor", "Length"] = "Factor",
        profile_resolution: InputInteger = 12,
        fill_caps: InputBoolean = False,
    ):
        super().__init__(
            **{
                "Curve": curve,
                "U Component": u_component,
                "Profile Resolution": profile_resolution,
                "Fill Caps": fill_caps,
            }
        )

    def _build_group(self, tree):
        curve = tree.inputs.geometry("Curve")
        u_component = tree.inputs.menu("U Component", optional_label=True)
        profile_resolution = tree.inputs.integer(
            "Profile Resolution", 12, min_value=3, max_value=512
        )
        fill_caps = tree.inputs.boolean("Fill Caps", False)
        mesh = tree.outputs.geometry("Mesh")
        uv_map = tree.outputs.vector("uv_map")

        _group = FallbackGeometry()
        menu_switch = g.MenuSwitch.integer(u_component, {"Factor": 0, "Length": 1})
        named_attribute = g.NamedAttribute.float("radius")
        _spline_length = g.SplineLength()
        spline_parameter = g.SplineParameter()
        index_switch = g.IndexSwitch.float(
            menu_switch.o.output, (spline_parameter.o.factor, spline_parameter.o.length)
        )
        capture = g.CaptureAttribute.point(
            geometry=g.CurveCircle(resolution=profile_resolution)
        )
        factor = capture.items.float("Factor", g.SplineParameter().o.factor)
        index = capture.items.integer("Index", g.SplineParameter().o.index)
        capture_1 = g.CaptureAttribute.point(geometry=curve)
        factor_1 = capture_1.items.float("Factor", index_switch)
        index_1 = capture_1.items.integer("Index", g.SplineParameter().o.index)
        length = capture_1.items.float("Length", g.SplineLength().o.length)
        switch = CheckEndFaceCorner(
            captured_index=index.output
        ).o.is_end_face_corner.switch.float(factor.output, 1.0)
        switch_1 = CheckEndFaceCorner(
            captured_index=index_1.output
        ).o.is_end_face_corner.switch.float(
            factor_1.output,
            g.IndexSwitch.float(menu_switch.o.output, (1.0, length.output)),
        )
        combine_xyz = g.CombineXYZ(x=switch_1, y=switch)
        (
            capture_1.o.geometry
            >> g.CurveToMesh(
                profile_curve=capture.o.geometry,
                scale=named_attribute.o.exists.switch.float(
                    1.0, named_attribute.o.attribute
                ),
                fill_caps=fill_caps,
            )
            >> mesh
        )

        combine_xyz >> uv_map

        u_component.default_value = "Factor"
