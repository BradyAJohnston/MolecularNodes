# Node group '.MN_utils_style_ribbon_nucleic' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MaterialSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMaterial,
    InputVector,
)
from ..atom_name import AtomName
from ..chain_id import ChainID
from ..curve_custom_profile import CurveCustomProfile
from ..curve_endpoint_values import CurveEndpointValues
from ..fallback_geometry import FallbackGeometry
from ..is_nucleic import IsNucleic
from ..offset_vector import OffsetVector
from ..unique_residue_id import UniqueResidueID
from .cleanup import Cleanup
from .curve_split_splines import CurveSplitSplines
from .curve_to_mesh_with_uvmap import CurveToMeshWithUVMap
from .mn_world_scale import MN_world_scale
from .sample_nucleic_base_values import SampleNucleicBaseValues
from .set_instancer import SetInstancer
from .smooth_by_angle import SmoothByAngle
from .vector_in_angstroms import VectorInAngstroms


class MN_utils_style_ribbon_nucleic(CustomGeometryGroup):
    """
    .MN_utils_style_ribbon_nucleic

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    material : InputMaterial
        Material to apply to the resulting geometry
    switch : InputBoolean
        Switch
    backbone_subdivisions : InputInteger
        Backbone Subdivisions
    backbone_resolution : InputInteger
        Backbone Resolution
    backbone_radius : InputFloat
        Backbone Radius
    backbone_shade_smooth : InputBoolean
        Backbone Shade Smooth
    backbone_scale : InputVector
        Backbone Scale
    base_geometry : InputGeometry
        Base Geometry
    base_scale : InputVector
        Base Scale
    base_resolution : InputInteger
        Base Resolution

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.material : MaterialSocket
        Material to apply to the resulting geometry
    i.switch : BooleanSocket
        Switch
    i.backbone_subdivisions : IntegerSocket
        Backbone Subdivisions
    i.backbone_resolution : IntegerSocket
        Backbone Resolution
    i.backbone_radius : FloatSocket
        Backbone Radius
    i.backbone_shade_smooth : BooleanSocket
        Backbone Shade Smooth
    i.backbone_scale : VectorSocket
        Backbone Scale
    i.base_geometry : GeometrySocket
        Base Geometry
    i.base_scale : VectorSocket
        Base Scale
    i.base_resolution : IntegerSocket
        Base Resolution

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    o.curve : GeometrySocket
        Curve
    """

    _name = ".MN_utils_style_ribbon_nucleic"
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry._mn_utils_style_ribbon_nucleic"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""
        switch: BooleanSocket
        """Switch"""
        backbone_subdivisions: IntegerSocket
        """Backbone Subdivisions"""
        backbone_resolution: IntegerSocket
        """Backbone Resolution"""
        backbone_radius: FloatSocket
        """Backbone Radius"""
        backbone_shade_smooth: BooleanSocket
        """Backbone Shade Smooth"""
        backbone_scale: VectorSocket
        """Backbone Scale"""
        base_geometry: GeometrySocket
        """Base Geometry"""
        base_scale: VectorSocket
        """Base Scale"""
        base_resolution: IntegerSocket
        """Base Resolution"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        curve: GeometrySocket
        """Curve"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        material: InputMaterial = None,
        switch: InputBoolean = False,
        backbone_subdivisions: InputInteger = 3,
        backbone_resolution: InputInteger = 8,
        backbone_radius: InputFloat = 2.0,
        backbone_shade_smooth: InputBoolean = True,
        backbone_scale: InputVector = None,
        base_geometry: InputGeometry = None,
        base_scale: InputVector = None,
        base_resolution: InputInteger = 6,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Material": material,
                "Switch": switch,
                "Backbone Subdivisions": backbone_subdivisions,
                "Backbone Resolution": backbone_resolution,
                "Backbone Radius": backbone_radius,
                "Backbone Shade Smooth": backbone_shade_smooth,
                "Backbone Scale": backbone_scale,
                "Base Geometry": base_geometry,
                "Base Scale": base_scale,
                "Base Resolution": base_resolution,
            }
        )

    def _build_group(self, tree):
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        material = tree.inputs.material(
            "Material", description="Material to apply to the resulting geometry"
        )
        switch = tree.inputs.boolean("Switch", False)
        with tree.inputs.panel("Backbone"):
            backbone_subdivisions = tree.inputs.integer(
                "Backbone Subdivisions", 3, min_value=1, max_value=10
            )
            backbone_resolution = tree.inputs.integer(
                "Backbone Resolution", 8, min_value=3, max_value=50
            )
            backbone_radius = tree.inputs.float(
                "Backbone Radius", 2.0, min_value=0.0, subtype="DISTANCE"
            )
            backbone_shade_smooth = tree.inputs.boolean("Backbone Shade Smooth", True)
            backbone_scale = tree.inputs.vector(
                "Backbone Scale", (1.0, 1.0, 1.0), subtype="XYZ"
            )
        with tree.inputs.panel("Base"):
            base_geometry = tree.inputs.geometry("Base Geometry")
            base_scale = tree.inputs.vector(
                "Base Scale", (2.5, 0.5, 7.0), min_value=-10_000.0, max_value=10_000.0
            )
            base_resolution = tree.inputs.integer(
                "Base Resolution", 6, min_value=3, max_value=512
            )
        geometry = tree.outputs.geometry("Geometry")
        curve = tree.outputs.geometry("Curve")

        with g.Frame("Slightly Extend Curve Ends"):
            endpoint_selection = g.EndpointSelection()
            group = CurveEndpointValues()
            group_1 = VectorInAngstroms(
                vector=OffsetVector(vector=g.CurveTangent(), offset=group.o.value),
                normalize=False,
                angstrom=group.o.value * -2.0,
            )
        _group_2 = CurveToMeshWithUVMap()
        remove_named_attribute = g.RemoveNamedAttribute(
            geometry=atoms, name="bond_type"
        )
        remove_named_attribute.node.warning_propagation = "ERRORS"
        capture = g.CaptureAttribute.point(geometry=remove_named_attribute)
        selection_1 = capture.items.boolean(
            "Selection", IsNucleic(and_=selection).o.selection
        )
        backbone_radius_1 = capture.items.float("Backbone Radius", backbone_radius)
        capture_1 = g.CaptureAttribute.point(
            geometry=(
                capture.o.geometry
                >> g.SeparateGeometry.point(selection=selection_1.output)
            ).o.selection
        )
        unique_group_id = capture_1.items.integer("Unique Group ID", UniqueResidueID())
        group_3 = SampleNucleicBaseValues(input=unique_group_id.output)
        capture_2 = g.CaptureAttribute.point(geometry=capture_1.o.geometry)
        base_valid = capture_2.items.boolean("base_valid", group_3.o.base_valid)
        capture_2.items.vector("base_pivot", group_3.o.base_pivot)
        base_z = capture_2.items.vector("base_Z", group_3.o.base_z)
        base_y = capture_2.items.vector("base_Y", group_3.o.base_y)
        base_position = capture_2.items.vector("base_position", group_3.o.base_position)
        base_color = capture_2.items.color("Base Color", group_3.o.base_color)
        with g.Frame("Delete between chains and distance too large"):
            points_to_curves = (
                capture_2.o.geometry
                >> g.MeshToPoints(
                    selection=g.Compare.integer.equal(AtomName(), 55),
                    position=base_position.output,
                    radius=0.05,
                )
                >> g.StoreNamedAttribute.point.integer(name="tmp_idx", value=g.Index())
                >> g.PointsToCurves(curve_group_id=ChainID())
            )
            group_4 = CurveSplitSplines(
                curve=points_to_curves,
                curve_normal="Minimum Twist",
                distance_cutoff=0.1,
            )
        set_curve_radius = group_4 >> g.SetCurveRadius(radius=backbone_radius_1.output)
        with g.Frame("Instance simple base cylinder"):
            axes_to_rotation = g.AxesToRotation(
                primary_axis=base_z.output, secondary_axis=base_y.output
            )
            store_named_attribute = (
                g.SetCurveNormal(
                    curve=set_curve_radius, normal=base_y.output, mode="Free"
                )
                >> g.StoreNamedAttribute.point.color(
                    name="vertex_color", value=base_color.output
                )
                >> g.StoreNamedAttribute.point.boolean(name="is_side_chain", value=True)
            )
            value = g.Value(1.0)
            cylinder = g.Cylinder(vertices=base_resolution, radius=value, depth=value)
            transform_geometry = (
                cylinder
                >> g.TransformGeometry(
                    translation=g.CombineXYZ(z=value.o.value / 2.0),
                    rotation=(0.0, 0.0, math.pi / 4),
                )
                >> g.StoreNamedAttribute.corner.vector(
                    name="uv_map", value=cylinder.o.uv_map
                )
                >> g.TransformGeometry(
                    scale=VectorInAngstroms(
                        vector=base_scale, normalize=False, angstrom=1.0
                    )
                )
            )
            instance_on_points = SetInstancer(
                geometry=store_named_attribute
            ) >> g.InstanceOnPoints(
                selection=base_valid.output,
                instance=FallbackGeometry(
                    geometry=base_geometry, fallback=transform_geometry
                ),
                rotation=axes_to_rotation,
            )
        store_named_attribute_1 = g.StoreNamedAttribute.point.boolean(
            g.SetCurveNormal(curve=set_curve_radius), name="is_backbone", value=True
        )
        set_position = g.SetPosition(
            geometry=store_named_attribute_1,
            selection=endpoint_selection,
            offset=group_1,
        )
        transform_geometry_1 = g.TransformGeometry(
            geometry=g.CurveCircle(resolution=4, radius=0.01),
            rotation=(math.pi / 2, math.pi / 4, -math.pi / 2),
        )
        remove_named_attribute_1 = (
            SmoothByAngle(
                mesh=instance_on_points, angle=5 * math.pi / 12, ignore_sharpness=True
            )
            >> g.SetMaterial(material=material)
            >> g.RemoveNamedAttribute(pattern_mode="Wildcard", name="tmp_*")
        )
        group_5 = CurveCustomProfile(
            curve=set_position,
            subdivisions=backbone_subdivisions,
            profile_type=switch.switch.menu("Custom Profile", "Default Profile"),
            socket_6=axes_to_rotation.o.rotation.rotate(
                (0.0, 0.0, -math.pi / 4), rotation_space="LOCAL"
            ),
            profile_scale=backbone_scale / MN_world_scale(),
            profile_curve=switch.switch.geometry(transform_geometry_1),
            profile_resolution=backbone_resolution,
            input_14=58.48539,
        )
        group_6 = Cleanup(
            geometry=SmoothByAngle(mesh=group_5, angle=math.pi / 3),
            color_source=store_named_attribute_1,
            material=material,
            shade_smooth=backbone_shade_smooth,
        )
        join_geometry = g.JoinGeometry(geometry=(group_6, remove_named_attribute_1))

        join_geometry >> geometry
        store_named_attribute_1 >> curve
