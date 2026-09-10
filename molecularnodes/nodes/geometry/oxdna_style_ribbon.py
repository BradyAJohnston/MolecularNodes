# Node-group asset 'oxDNA Style Ribbon' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ColorSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
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
    InputInteger,
    InputMaterial,
    InputVector,
)
from ._shared.mn_units import MNUnits
from ._shared.set_instancer import SetInstancer
from .angstrom_to_world import AngstromToWorld
from .color_res_name import ColorResName
from .oxdna_normal import OxDNANormal
from .oxdna_offset import OxDNAOffset
from .oxdna_rotation import OxDNARotation
from .set_color import SetColor


class Utils_oxdna_base(CustomGeometryGroup):
    _name = ".utils_oxdna_base"
    _tree_properties = {"node_tool_idname": "geometry._utils_oxdna_base"}

    def _build_group(self, tree):
        value = tree.inputs.float("Value", 0.5, min_value=-10_000.0, max_value=10_000.0)
        value_1 = tree.inputs.float(
            "Value", 0.5, min_value=-10_000.0, max_value=10_000.0
        )
        value_2 = tree.inputs.float(
            "Value", 0.5, min_value=-10_000.0, max_value=10_000.0
        )
        geometry = tree.outputs.geometry("Geometry")

        group = MNUnits(value=value_1)
        (
            g.Cylinder(
                radius=MNUnits(value=value).o.angstrom,
                depth=group.o.angstrom,
                vertices=4,
            )
            >> g.TransformGeometry(translation=g.CombineXYZ(z=group.o.angstrom / 2.0))
            >> g.TransformGeometry(rotation=(0.0, 0.0, math.pi / 4))
            >> g.TransformGeometry(
                scale=g.CombineXYZ(x=MNUnits(value=value_2).o.angstrom, y=1.0, z=1.0)
            )
            >> geometry
        )


class OxDNAStyleRibbon(AssetGeometryGroup):
    """
    oxDNA Style Ribbon

    Parameters
    ----------
    atoms : InputGeometry
        Atoms
    selection : InputBoolean
        Selection of atoms to apply this node to
    backbone_resolution : InputInteger
        Backbone Resolution
    backbone_subdivisions : InputInteger
        Backbone Subdivisions
    backbone_radius : InputFloat
        Backbone Radius
    a : InputColor
        A
    c : InputColor
        C
    g : InputColor
        G
    t_u : InputColor
        T / U
    base_scale : InputVector
        Base Scale
    shade_smooth : InputBoolean
        Shade Smooth
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.atoms : GeometrySocket
        Atoms
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.backbone_resolution : IntegerSocket
        Backbone Resolution
    i.backbone_subdivisions : IntegerSocket
        Backbone Subdivisions
    i.backbone_radius : FloatSocket
        Backbone Radius
    i.a : ColorSocket
        A
    i.c : ColorSocket
        C
    i.g : ColorSocket
        G
    i.t_u : ColorSocket
        T / U
    i.base_scale : VectorSocket
        Base Scale
    i.shade_smooth : BooleanSocket
        Shade Smooth
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "oxDNA Style Ribbon"
    _asset_name = "oxDNA Style Ribbon"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.mn_oxdna_style_ribbon"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atoms"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        backbone_resolution: IntegerSocket
        """Backbone Resolution"""
        backbone_subdivisions: IntegerSocket
        """Backbone Subdivisions"""
        backbone_radius: FloatSocket
        """Backbone Radius"""
        a: ColorSocket
        """A"""
        c: ColorSocket
        """C"""
        g: ColorSocket
        """G"""
        t_u: ColorSocket
        """T / U"""
        base_scale: VectorSocket
        """Base Scale"""
        shade_smooth: BooleanSocket
        """Shade Smooth"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""

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
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        backbone_resolution: InputInteger = 6,
        backbone_subdivisions: InputInteger = 1,
        backbone_radius: InputFloat = 2.0,
        a: InputColor = None,
        c: InputColor = None,
        g: InputColor = None,
        t_u: InputColor = None,
        base_scale: InputVector = None,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Backbone Resolution": backbone_resolution,
                "Backbone Subdivisions": backbone_subdivisions,
                "Backbone Radius": backbone_radius,
                "A": a,
                "C": c,
                "G": g,
                "T / U": t_u,
                "Base Scale": base_scale,
                "Shade Smooth": shade_smooth,
                "Material": material,
            }
        )

    def _build_group(self, tree):
        atoms = tree.inputs.geometry("Atoms")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        with tree.inputs.panel("Backbone", default_closed=True):
            backbone_resolution = tree.inputs.integer(
                "Backbone Resolution", 6, min_value=3, max_value=512
            )
            backbone_subdivisions = tree.inputs.integer(
                "Backbone Subdivisions", 1, min_value=1
            )
            backbone_radius = tree.inputs.float(
                "Backbone Radius", 2.0, min_value=0.0, max_value=10_000.0
            )
        with tree.inputs.panel("Base", default_closed=True):
            a = tree.inputs.color("A", (0.2746774, 0.5457247, 0.799103, 1.0))
            c_ = tree.inputs.color("C", (0.294582, 0.8, 0.187789, 1.0))
            g_ = tree.inputs.color("G", (0.8, 0.236614, 0.167417, 1.0))
            t_u = tree.inputs.color("T / U", (0.8, 0.269803, 0.526898, 1.0))
            base_scale = tree.inputs.vector(
                "Base Scale", (1.0, 1.0, 1.0), subtype="XYZ"
            )
        with tree.inputs.panel("Material", default_closed=True):
            shade_smooth = tree.inputs.boolean("Shade Smooth", True)
            material = tree.inputs.material(
                "Material", description="Material to apply to the resulting geometry"
            )
        geometry = tree.outputs.geometry("Geometry")

        capture = g.CaptureAttribute.point(
            geometry=g.SeparateGeometry.point(atoms, selection).o.selection
        )
        value = capture.items.vector("Value", OxDNAOffset())
        with g.Frame("Colored bases"):
            group = SetColor(
                atoms=capture.o.geometry,
                color=ColorResName(a=a, c=c_, g=g_, t=t_u, ra=a, rc=c_, rg=g_, ru=t_u),
            )
            instance_on_points = SetInstancer(geometry=group) >> g.InstanceOnPoints(
                instance=Utils_oxdna_base(
                    _named_links=[("Value", 4.139999), ("Value", 5.42), ("Value", 3.32)]
                ),
                rotation=OxDNARotation(),
                scale=base_scale,
            )
        set_position = g.SetPosition(geometry=capture.o.geometry, offset=value.output)
        with g.Frame("Backbone ribbon"):
            set_spline_type = (
                g.MeshToCurve(mesh=set_position)
                >> g.SetCurveRadius(radius=AngstromToWorld(angstrom=backbone_radius))
                >> g.SetSplineType.bezier()
            )
            set_shade_smooth = (
                g.SetHandleType(curve=set_spline_type)
                >> g.SetSplineResolution(resolution=backbone_subdivisions)
                >> g.CurveToMesh(
                    profile_curve=g.CurveCircle(resolution=backbone_resolution),
                    scale=g.Radius(),
                    fill_caps=True,
                )
                >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
            )
        group_1 = SetInstancer(geometry=set_position)
        with g.Frame("Base stem"):
            group_2 = Utils_oxdna_base(
                _named_links=[("Value", 1.4999993), ("Value", 6.119999), ("Value", 5.0)]
            )
            instance_on_points_1 = group_1 >> g.InstanceOnPoints(
                instance=group_2,
                rotation=g.AxesToRotation(
                    primary_axis=value.output * -1.0, secondary_axis=OxDNANormal()
                ),
            )
        (
            g.JoinGeometry(
                geometry=(set_shade_smooth, instance_on_points_1, instance_on_points)
            )
            >> g.SetMaterial(material=material)
            >> geometry
        )
        _group_3 = OxDNARotation()


ASSET = OxDNAStyleRibbon

ASSET_METADATA = {
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}
