# Node-group asset "Style Cartoon" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MaterialSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMaterial,
    InputMenu,
    InputVector,
)
from ._shared.cleanup import Cleanup
from ._shared.curve_split_splines import CurveSplitSplines
from ._shared.mn_units import MNUnits
from ._shared.mn_utils_style_ribbon_nucleic import MN_utils_style_ribbon_nucleic
from ._shared.vector_in_angstroms import VectorInAngstroms
from .atoms_to_ca_curves import AtomsToCACurves
from .boolean_run_trim import BooleanRunTrim
from .curve_custom_profile import CurveCustomProfile
from .curve_endpoint_values import CurveEndpointValues
from .curve_offset_dot import CurveOffsetDot
from .curve_rotation import CurveRotation
from .dihedral_phi import DihedralPhi
from .dihedral_psi import DihedralPsi
from .evaluate_on_atoms import EvaluateOnAtoms
from .expand_boolean import ExpandBoolean
from .is_helix import IsHelix
from .is_loop import IsLoop
from .is_sheet import IsSheet
from .offset_boolean import OffsetBoolean
from .offset_color_attribute import OffsetColorAttribute
from .offset_point_along_curve import OffsetPointAlongCurve
from .offset_rotation import OffsetRotation
from .offset_vector import OffsetVector
from .separate_polymers import SeparatePolymers
from .set_color import SetColor
from .sub_group_info import SubGroupInfo


class Tmp_ss_attributes(CustomGeometryGroup):
    _name = ".tmp_ss_attributes"
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry._tmp_ss_attributes"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        tmp_ss_is_first = tree.outputs.boolean("tmp_ss_is_first")
        tmp_ss_is_last = tree.outputs.boolean("tmp_ss_is_last")
        tmp_ss_size = tree.outputs.integer("tmp_ss_size")
        tmp_idx = tree.outputs.integer("tmp_idx")
        tmp_idx_curve = tree.outputs.integer("tmp_idx_curve")
        tmp_curve_normal = tree.outputs.vector("tmp_curve_normal")
        tmp_curve_tangent = tree.outputs.vector("tmp_curve_tangent")

        named_attribute = g.NamedAttribute.boolean("tmp_ss_is_first")
        named_attribute_1 = g.NamedAttribute.boolean("tmp_ss_is_last")
        named_attribute_2 = g.NamedAttribute.integer("tmp_ss_size")
        named_attribute_3 = g.NamedAttribute.integer("tmp_idx")
        named_attribute_4 = g.NamedAttribute.integer("tmp_idx_curve")
        named_attribute_5 = g.NamedAttribute.vector("tmp_curve_normal")
        named_attribute_6 = g.NamedAttribute.vector("tmp_curve_tangent")

        named_attribute >> tmp_ss_is_first
        named_attribute_1 >> tmp_ss_is_last
        named_attribute_2 >> tmp_ss_size
        named_attribute_3 >> tmp_idx
        named_attribute_4 >> tmp_idx_curve
        named_attribute_5 >> tmp_curve_normal
        named_attribute_6 >> tmp_curve_tangent


class SampleFromCACurve(CustomGeometryGroup):
    _name = ".Sample from CA curve"
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry._sample_from_ca_curve"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        ca_curve = tree.inputs.geometry("CA Curve")
        offset = tree.inputs.float(
            "Offset", 0.0, min_value=-10_000.0, max_value=10_000.0
        )
        position = tree.outputs.vector("Position")
        tangent = tree.outputs.vector("Tangent")
        normal = tree.outputs.vector("Normal")

        group = Tmp_ss_attributes()
        sample_curve = g.SampleCurve(
            curves=ca_curve,
            factor=OffsetPointAlongCurve(
                point_index=group.o.tmp_idx, offset=offset
            ).o.factor,
            curve_index=group.o.tmp_idx_curve,
        )
        capture = g.CaptureAttribute.point(geometry=ca_curve)
        position_1 = capture.items.vector("Position", sample_curve.o.position)
        tangent_1 = capture.items.vector("Tangent", sample_curve.o.tangent)
        normal_1 = capture.items.vector("Normal", sample_curve.o.normal)
        sample_index = g.SampleIndex(
            geometry=capture.o.geometry,
            value=position_1.output,
            index=group.o.tmp_idx,
            data_type="FLOAT_VECTOR",
        )
        sample_index_1 = g.SampleIndex(
            geometry=capture.o.geometry,
            value=normal_1.output,
            index=group.o.tmp_idx,
            data_type="FLOAT_VECTOR",
        )
        sample_index_2 = g.SampleIndex(
            geometry=capture.o.geometry,
            value=tangent_1.output,
            index=group.o.tmp_idx,
            data_type="FLOAT_VECTOR",
        )

        sample_index >> position
        sample_index_2 >> tangent
        sample_index_1 >> normal


class FixLoopAlignmentIntoAH(CustomGeometryGroup):
    _name = ".Fix Loop Alignment into AH"
    _tree_properties = {"node_tool_idname": "geometry._fix_loop_alignment_into_ah"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        curve = tree.inputs.geometry("Curve")
        tangent = tree.inputs.vector("Tangent", (0.0, 0.0, 0.0), hide_value=True)
        geometry = tree.outputs.geometry("Geometry")

        capture = g.CaptureAttribute.point(geometry=curve)
        selection = capture.items.boolean("Selection", IsHelix().o.selection)
        with g.Frame("Subdivide the first and last segments going into AH"):
            boolean_math = g.BooleanMath.subtract(
                g.EndpointSelection(start_size=0, end_size=2),
                g.EndpointSelection(start_size=0),
            )
            boolean_math_1 = ExpandBoolean(
                boolean=selection.output, expand=1
            ).o.boolean & (boolean_math.o.boolean | g.EndpointSelection(end_size=0))
            subdivide_curve = capture.o.geometry >> g.SubdivideCurve(
                cuts=boolean_math_1.switch.integer(true=1)
            )
        with g.Frame("Selection for second point and second last point"):
            boolean_math_2 = g.BooleanMath.subtract(
                g.EndpointSelection(start_size=0, end_size=2),
                g.EndpointSelection(start_size=0),
            )
            boolean_math_3 = g.BooleanMath.subtract(
                g.EndpointSelection(start_size=2, end_size=0),
                g.EndpointSelection(end_size=0),
            )
            boolean_math_4 = (
                boolean_math_2.o.boolean | boolean_math_3
            ) & selection.output
        capture_1 = g.CaptureAttribute.point(geometry=subdivide_curve)
        value = capture_1.items.integer(
            "Value",
            CurveEndpointValues(
                start_size=2, start_value=-1, end_size=2, end_value=1
            ).o.value,
        )
        group = VectorInAngstroms(
            vector=OffsetVector(vector=tangent, offset=value.output),
            normalize=False,
            angstrom=value.output * -0.5,
        )
        set_position = capture_1.o.geometry >> g.SetPosition(
            selection=boolean_math_4,
            position=OffsetVector(vector=g.Position(), offset=value.output),
            offset=group,
        )
        with g.Frame("Move loop ends slightly inside of SS"):
            (
                set_position
                >> g.SetPosition(
                    selection=g.EndpointSelection(),
                    offset=VectorInAngstroms(
                        vector=tangent, normalize=False, angstrom=value.output * 0.3
                    ),
                )
                >> geometry
            )


class CAToLoops(CustomGeometryGroup):
    _name = ".CA to loops"
    _tree_properties = {"node_tool_idname": "geometry._ca_to_loops"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        subdivisions = tree.inputs.integer("Subdivisions", 12, min_value=1)
        radius = tree.inputs.float(
            "Radius", 0.0, min_value=-10_000.0, max_value=10_000.0
        )
        as_cylinders = tree.inputs.boolean("As Cylinders", False)
        profile_resolution = tree.inputs.integer(
            "Profile Resolution", 6, min_value=3, max_value=512
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        endpoint_selection = g.EndpointSelection(end_size=0)
        group = IsLoop()
        with g.Frame("Find where it transitions direction from one SS to another"):
            group_1 = IsHelix()
            group_2 = IsSheet()
            boolean_math = ExpandBoolean(
                boolean=group_1.o.selection, expand=1
            ).o.boolean & group_2.o.selection | group_1.o.selection & ExpandBoolean(
                boolean=group_2.o.selection, expand=1
            )
        group_3 = IsHelix(and_=g.EndpointSelection().o.selection & as_cylinders)
        group_4 = Tmp_ss_attributes()
        group_5 = SampleFromCACurve(**{"CA Curve": geometry}, Offset=-0.1)
        group_6 = SampleFromCACurve(**{"CA Curve": geometry}, Offset=0.1)
        with g.Frame("Expand selection by 1 so the ribbon stops where the SS starts"):
            boolean_math_1 = (
                ExpandBoolean(boolean=group.o.selection, expand=1).o.boolean
                | boolean_math
            )
        with g.Frame("catch direct change from one to another (needs improving)"):
            _boolean_math_2 = IsHelix().o.selection & OffsetBoolean(
                boolean=IsSheet().o.selection, offset=1
            ) | IsSheet().o.selection & OffsetBoolean(
                boolean=IsHelix().o.selection, offset=-1
            )
        capture = g.CaptureAttribute.point(geometry=geometry)
        boolean = capture.items.boolean("Boolean", boolean_math_1)
        subdivisions_1 = capture.items.integer("Subdivisions", subdivisions)
        radius_1 = capture.items.float("Radius", radius)
        group_7 = VectorInAngstroms(
            vector=g.Normal(legacy_corner_normals=True).o.normal,
            normalize=False,
            angstrom=2.0,
        )
        rotate_rotation = CurveRotation().o.rotation.rotate(
            (0.0, 0.0, 0.0), rotation_space="LOCAL"
        )
        capture_1 = g.CaptureAttribute.point(geometry=capture.o.geometry)
        rotation = capture_1.items.rotation("Rotation", rotate_rotation)
        with g.Frame("Don't resample when directly from one SS to another"):
            compare = g.SplineLength().o.point_count > 2
        switch = group_4.o.tmp_ss_is_last.switch.float(
            group_4.o.tmp_ss_is_first.switch.float(true=0.04), -0.09
        )
        group_8 = CurveSplitSplines(
            curve=capture_1.o.geometry,
            selection=boolean.output,
            distance_cutoff=0.05,
            rotation=rotation.output,
            offset_amount=switch,
        )
        switch_1 = endpoint_selection.o.selection.switch.vector(
            group_6.o.position, group_5.o.position
        )
        switch_2 = endpoint_selection.o.selection.switch.vector(
            SampleFromCACurve(**{"CA Curve": geometry}, Offset=0.2).o.tangent,
            SampleFromCACurve(**{"CA Curve": geometry}, Offset=-0.2).o.tangent,
        )
        capture_2 = g.CaptureAttribute.point(geometry=group_8)
        capture_2.items.vector("Position", switch_1)
        tangent = capture_2.items.vector("Tangent", switch_2)
        set_position = (
            capture_2.o.geometry
            >> g.SetCurveNormal(
                selection=group_3.o.selection,
                normal=endpoint_selection.o.selection.switch.vector(
                    group_6.o.normal, group_5.o.normal
                ),
                mode="Free",
            )
            >> g.SetPosition(
                selection=group_3.o.selection, position=switch_1, offset=group_7
            )
        )
        group_9 = SetColor(
            atoms=set_position,
            selection=compare,
            color=OffsetColorAttribute(offset=CurveEndpointValues().o.value),
        )
        set_curve_radius = (
            g.SetHandleType(
                curve=FixLoopAlignmentIntoAH(Curve=group_9, Tangent=tangent.output)
            )
            >> g.SetCurveNormal()
            >> g.SetCurveRadius(radius=radius_1.output)
        )
        (
            CurveCustomProfile(
                curve=set_curve_radius,
                subdivisions=subdivisions_1.output,
                profile_type="Default Profile",
                socket_6=rotate_rotation,
                profile_resolution=profile_resolution,
                input_14=0.0,
            )
            >> g.StoreNamedAttribute.point.integer(name="sec_struct", value=3)
            >> geometry_1
        )


class SplitCurves(CustomGeometryGroup):
    _name = "Split Curves"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        curves = tree.inputs.geometry("Curves")
        selection = tree.inputs.boolean(
            "Selection", False, description="The calculated selection"
        )
        resolution = tree.inputs.integer("Resolution", 12, min_value=1)
        expand = tree.inputs.integer("Expand", 1, min_value=-2147483647)
        curve = tree.outputs.geometry("Curve")

        _group = DihedralPhi()
        compare = DihedralPsi().o.psi > -3.080002
        with g.Frame("Switch back and forth for beta-sheet fixing"):
            boolean_math = IsSheet().o.selection & (
                g.Index().o.index - g.AccumulateField.point.integer(compare).o.leading
            ).modulo(2)
            _switch = boolean_math.switch.float(-1.0, 1.0)
        _group_1 = DihedralPsi()
        accumulate_field = g.AccumulateField.point.integer(
            g.BooleanMath.subtract(
                IsSheet().o.selection, DihedralPhi().o.phi > 2.719998
            )
        )
        curve_to_points = curves >> g.SetSplineType() >> g.CurveToPoints.evaluated()
        switch_1 = g.Compare.integer.not_equal(expand, 0).o.result.switch.boolean(
            selection, ExpandBoolean(boolean=selection, expand=expand)
        )
        capture = g.CaptureAttribute.point(geometry=curve_to_points)
        selection_1 = capture.items.boolean("Selection", switch_1)
        vector_math = curve_to_points.o.normal * g.Switch.float(
            accumulate_field.o.leading.modulo(2), -1.0, 1.0
        )
        accumulate_field_1 = g.AccumulateField.point.integer(
            g.Position().o.position.distance(OffsetVector(offset=-1)) > 5.0
        )
        capture_1 = g.CaptureAttribute.point(geometry=capture.o.geometry)
        group_id = capture_1.items.integer(
            "Group ID",
            SubGroupInfo(
                sub_group_id=selection_1.output, group_id=accumulate_field_1.o.trailing
            ).o.group_id,
        )
        set_spline_type = (
            g.SeparateGeometry.point(capture_1.o.geometry, selection_1.output)
            >> g.PointsToCurves(curve_group_id=group_id.output)
            >> g.SeparateGeometry.spline(selection=g.SplineLength().o.point_count > 1)
            >> g.SetCurveNormal(normal=vector_math, mode="Free")
            >> g.SetSplineResolution(resolution=resolution)
            >> g.SetSplineType.bezier()
        )
        g.SetHandleType(curve=set_spline_type) >> curve
        viewer = g.Viewer()
        capture_1.o.geometry >> viewer
        selection >> viewer


class NodeGroup(CustomGeometryGroup):
    _name = "NodeGroup"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        curves = tree.inputs.geometry(
            "Curves",
            description="Geometry to evaluate the given fields and store the resulting attributes on. All geometry types except volumes are supported",
        )
        _reset_curve_normal = tree.inputs.boolean(
            "Reset Curve Normal", False, structure_type="SINGLE", force_non_field=True
        )
        rotation = tree.inputs.rotation("Rotation", (0.0, 0.0, 0.0))
        menu = tree.inputs.menu("Menu", optional_label=True)
        resolution = tree.inputs.integer(
            "Resolution",
            12,
            description="Number of points on the circle",
            min_value=3,
            max_value=512,
        )
        radius = tree.inputs.float(
            "Radius",
            0.1,
            description="Distance of the points from the origin",
            min_value=0.0,
            subtype="DISTANCE",
        )
        scale = tree.inputs.vector("Scale", (1.0, 2.5, 1.0), subtype="XYZ")
        geometry = tree.outputs.geometry("Geometry")

        mix = g.Mix(b_float=1.0, clamp_factor=True)
        multiply_matrices = g.MultiplyMatrices(
            matrix=g.CombineTransform(rotation=rotation, scale=scale),
            matrix_001=g.CombineTransform(rotation=(0.0, 0.0, math.pi / 4)),
        )
        capture = g.CaptureAttribute.point(geometry=curves)
        rotation_1 = capture.items.rotation("Rotation", CurveRotation())
        capture.items.vector("Position", g.Position())
        sample_curve = g.SampleCurve(
            curves=capture.o.geometry,
            value=rotation_1.output,
            length=g.SplineParameter().o.length + 0.0,
            curve_index=g.CurveOfPoint().o.curve_index,
            mode="LENGTH",
            data_type="QUATERNION",
        )
        capture_1 = g.CaptureAttribute.point(geometry=capture.o.geometry)
        rotation_2 = capture_1.items.rotation("Rotation", sample_curve.o.value)
        position = capture_1.items.vector("Position", sample_curve.o.position)
        mix_1 = g.Mix(
            a_rotation=rotation_2.output,
            b_rotation=OffsetRotation(rotation=rotation_2.output, offset=1),
            factor_float=0.0,
            data_type="ROTATION",
            clamp_factor=True,
        )
        switch = g.EndpointSelection(end_size=0).o.selection.switch.rotation(
            rotation_2.output, mix_1.o.result_rotation
        )
        switch_1 = IsSheet().o.selection.switch.rotation(
            switch, switch.rotate((0.0, 0.0, math.pi / 4), rotation_space="LOCAL")
        )
        blur_attribute = g.BlurAttribute.vector(
            g.RotateVector(rotation=switch_1, vector=(1.0, 0.0, 0.0)),
            84,
            g.EndpointSelection(start_size=4, end_size=3),
        )
        set_curve_normal = (
            capture_1.o.geometry
            >> g.SetPosition(position=position.output)
            >> g.SetCurveNormal(normal=blur_attribute, mode="Free")
        )
        set_curve_normal.node.mute = True
        transform_geometry = g.TransformGeometry(
            geometry=g.CurveCircle(resolution=4, radius=0.1),
            transform=multiply_matrices,
            mode="Matrix",
            rotation=(0.0, 0.0, math.pi / 4),
        )
        fillet_curve = g.SetHandleType(
            curve=g.SetSplineType.bezier(transform_geometry)
        ) >> g.FilletCurve(
            radius=g.Mix(
                factor_float=mix.o.result_float,
                a_float=0.001,
                b_float=0.05,
                clamp_factor=True,
            ).o.result_float,
            limit_radius=True,
        )
        menu_switch = g.MenuSwitch.geometry(
            menu,
            {
                "Sharp": transform_geometry,
                "Smooth": (mix.o.result_float > 0.0).switch.geometry(
                    transform_geometry, fillet_curve
                ),
                "Round": g.CurveCircle(resolution=resolution, radius=radius),
            },
        )
        triangulate = (
            set_curve_normal
            >> g.CurveToMesh(profile_curve=menu_switch, fill_caps=True)
            >> g.SetShadeSmooth.face(shade_smooth=False)
            >> g.Triangulate()
        )
        triangulate.node.mute = True

        triangulate >> geometry

        menu.default_value = "Smooth"


class TweakArrowHeads(CustomGeometryGroup):
    _name = ".Tweak Arrow Heads"
    _tree_properties = {"node_tool_idname": "geometry._tweak_arrow_heads"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        rounded = tree.inputs.boolean("Rounded", False)
        input = tree.inputs.float("Input", 0.0)
        output = tree.outputs.geometry("Output")

        group = VectorInAngstroms(
            vector=Tmp_ss_attributes().o.tmp_curve_tangent,
            normalize=False,
            angstrom=input.map_range(to_max=-0.2),
        )
        vector_math = g.Normal(legacy_corner_normals=True).o.normal.dot(
            Tmp_ss_attributes().o.tmp_curve_normal
        )
        boolean_math = (input < 1.0) & (
            g.EvaluateOnDomain.face.float(
                g.Compare.float.equal(abs(vector_math), 0.0, 0.5)
            )
            > 0.44
        )
        capture = g.CaptureAttribute.face(geometry=geometry)
        boolean = capture.items.boolean("Boolean", boolean_math)
        normal = capture.items.vector(
            "Normal", g.Normal(legacy_corner_normals=True).o.normal
        )
        extrude_mesh = capture.o.geometry >> g.ExtrudeMesh(
            selection=boolean.output,
            offset=normal.output,
            offset_scale=input.map_range(to_max=MNUnits(value=0.8900002).o.angstrom),
            individual=False,
        )
        group_1 = VectorInAngstroms(
            vector=normal.output,
            normalize=False,
            angstrom=(input > 0.9).switch.float(
                input.map_range(from_max=3.77, to_max=-0.14)
            ),
        )
        set_position = g.SetPosition(
            geometry=g.SetPosition(
                geometry=extrude_mesh, selection=extrude_mesh.o.top, offset=group_1
            ),
            selection=extrude_mesh.o.top,
            offset=group,
        )
        rounded.switch.geometry(set_position, extrude_mesh) >> output


class CAToSheet(CustomGeometryGroup):
    _name = ".CA to sheet"
    _tree_properties = {"node_tool_idname": "geometry._ca_to_sheet"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        curve = tree.inputs.geometry("Curve")
        profile_resolution = tree.inputs.integer(
            "Profile Resolution", 4, min_value=3, max_value=512
        )
        thickness = tree.inputs.float(
            "Thickness", 2.34, min_value=-10_000.0, max_value=10_000.0
        )
        width = tree.inputs.float(
            "Width", 0.56, min_value=-10_000.0, max_value=10_000.0
        )
        subdivisions = tree.inputs.integer("Subdivisions", 6, min_value=1)
        rounded = tree.inputs.boolean("Rounded", False)
        arrows = tree.inputs.boolean("Arrows", False)
        geometry = tree.outputs.geometry("Geometry")

        with g.Frame("alternating flipping of normal to smoothen out sheets"):
            math_1 = g.Math(
                value=g.Index(),
                value_001=CurveOffsetDot(
                    offset=1, threshold_direction="Greater Than", threshold_cutoff=-0.3
                ).o.leading,
            )
            axis_angle_to_rotation = g.AxisAngleToRotation(
                axis=g.CurveTangent(),
                angle=math_1.o.value.wrap(0.0, 2.0) * 3.14159 + g.Math.to_radians(30.0),
            )
            blur_attribute = g.BlurAttribute.vector(
                g.Normal(legacy_corner_normals=True).o.normal.rotate(
                    axis_angle_to_rotation
                )
            )
        capture = g.CaptureAttribute.point(
            geometry=g.SetCurveNormal(curve=curve, normal=blur_attribute, mode="Free")
        )
        resolution = capture.items.integer("Resolution", subdivisions)
        thickness_1 = capture.items.float("Thickness", thickness)
        width_1 = capture.items.float("Width", width)
        vector_math = g.CurveTangent().o.tangent * (
            CurveEndpointValues(start_value=-1, end_value=1).o.value
            * MNUnits(value=0.2).o.angstrom
        )
        set_spline_type = (
            CurveSplitSplines(
                curve=capture.o.geometry,
                selection=BooleanRunTrim(boolean=IsSheet().o.selection, size=2),
                distance_cutoff=0.05,
                offset_spline_type="Poly",
            )
            >> g.SetPosition(offset=vector_math)
            >> g.SetSplineType.bezier()
        )
        capture_1 = g.CaptureAttribute.point(
            geometry=g.SetHandleType(curve=set_spline_type)
        )
        arrow_mask = capture_1.items.float(
            "Arrow Mask",
            (g.SplineLength().o.point_count - 1.0).max(1.0)
            - g.SplineParameter().o.index,
        )
        with g.Frame("Adjustment for arrowheads"):
            switch = arrows.switch.float(
                1.0,
                arrow_mask.output.map_range(
                    to_min=width_1.output * -1.0, to_max=-0.16000003
                ),
            )
        store_named_attribute = (
            capture_1.o.geometry
            >> g.StoreNamedAttribute.point.vector(
                name="tmp_curve_normal",
                value=g.Normal(legacy_corner_normals=True).o.normal,
            )
            >> g.StoreNamedAttribute.point.vector(
                name="tmp_curve_tangent", value=g.CurveTangent()
            )
        )
        combine_xyz = g.CombineXYZ(
            x=thickness_1.output,
            y=rounded.switch.float(width_1.output + switch, width_1.output),
        )
        group = CurveCustomProfile(
            curve=store_named_attribute,
            subdivisions=resolution.output,
            socket_6=CurveRotation(),
            profile_scale=combine_xyz,
            profile_resolution=profile_resolution,
            input_14=0.0,
        )
        (
            g.BooleanMath.subtract(arrows, rounded).o.boolean.switch.geometry(
                group,
                TweakArrowHeads(
                    Geometry=group, Rounded=rounded, Input=arrow_mask.output
                ),
            )
            >> geometry
        )


class BooleanShrink(CustomGeometryGroup):
    _name = "Boolean Shrink"
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.boolean_shrink"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        boolean = tree.inputs.boolean("Boolean", False, hide_value=True)
        shrink = tree.inputs.integer("Shrink", 0, min_value=-2147483647)
        boolean_1 = tree.outputs.boolean("Boolean")

        boolean_math = OffsetBoolean(
            boolean=boolean, offset=shrink
        ).o.boolean & OffsetBoolean(boolean=boolean, offset=shrink * -1.0)
        (boolean & boolean_math) >> boolean_1


class CAToHelix(CustomGeometryGroup):
    _name = ".CA to helix"
    _tree_properties = {"node_tool_idname": "geometry._ca_to_helix"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        curve = tree.inputs.geometry("Curve")
        boolean = tree.inputs.boolean("Boolean", False)
        thickness = tree.inputs.float(
            "Thickness", 2.34, min_value=-10_000.0, max_value=10_000.0
        )
        width = tree.inputs.float(
            "Width", 0.56, min_value=-10_000.0, max_value=10_000.0
        )
        profile_resolution = tree.inputs.integer(
            "Profile Resolution", 4, min_value=3, max_value=512
        )
        subdivisions = tree.inputs.integer("Subdivisions", 6, min_value=1)
        geometry = tree.outputs.geometry("Geometry")

        with g.Frame("Creating Alpha-helix Geometry"):
            _group = CurveRotation()
            rotate_rotation = (
                CurveRotation()
                .o.rotation.rotate((math.pi / 9, 0.0, 0.0), rotation_space="LOCAL")
                .rotate(
                    CurveOffsetDot(offset=1, threshold_cutoff=-0.7).o.rotation,
                    rotation_space="LOCAL",
                )
            )
            align_rotation_to_vector = g.AlignRotationToVector(
                rotation=OffsetRotation(
                    rotation=rotate_rotation, offset=CurveEndpointValues().o.value
                ),
                vector=g.CurveTangent(),
                pivot_axis="Y",
            )
            sample_curve = g.SampleCurve(
                curves=curve,
                value=CurveRotation(),
                length=g.SplineParameter().o.length + 0.2,
                curve_index=g.CurveOfPoint().o.curve_index,
                mode="LENGTH",
                data_type="QUATERNION",
            )
            capture = g.CaptureAttribute.point(geometry=curve)
            capture.node.mute = True
            value = capture.items.rotation("Value", sample_curve.o.value)
            position = capture.items.vector("Position", sample_curve.o.position)
            set_position = capture.o.geometry >> g.SetPosition(position=position.output)
            set_position.node.mute = True
            set_spline_type = (
                CurveSplitSplines(
                    curve=set_position,
                    selection=IsHelix(and_=~boolean).o.selection,
                    distance_cutoff=0.05,
                    rotation=value.output,
                    offset_amount=0.15,
                )
                >> g.SetSplineType.bezier()
            )
            group_1 = CurveCustomProfile(
                curve=g.SetHandleType(curve=set_spline_type),
                subdivisions=subdivisions,
                socket_6=align_rotation_to_vector,
                profile_scale=g.CombineXYZ(x=thickness, y=width),
                profile_resolution=profile_resolution,
                input_14=0.0,
            )
        with g.Frame("Creating Helix Cylinders"):
            boolean_math = ~g.EndpointSelection().o.selection
            _group_2 = ExpandBoolean()
            _group_3 = CurveCustomProfile(
                profile_scale=(2.65, 2.65, 2.65), profile_resolution=8, input_14=0.0
            )
            group_4 = IsHelix(and_=boolean)
            named_attribute = g.NamedAttribute.float("radius")
            group_5 = CurveSplitSplines(
                curve=g.SetCurveRadius(curve=curve, radius=0.085),
                selection=group_4.o.selection,
                distance_cutoff=0.06,
            )
            set_curve_radius = (
                g.SetPosition(
                    geometry=group_5,
                    offset=VectorInAngstroms(
                        vector=g.Normal(legacy_corner_normals=True).o.normal,
                        angstrom=2.4,
                    ),
                )
                >> g.SetSplineType()
                >> g.SetPosition(
                    selection=boolean_math,
                    position=g.BlurAttribute.vector(g.Position(), 2, boolean_math),
                )
                >> g.SetCurveNormal()
                >> g.ResampleCurve(length=MNUnits(value=2.0).o.angstrom, mode="Length")
                >> g.SetSplineType.bezier()
                >> g.SetCurveRadius(radius=MNUnits(value=width * 4.0).o.angstrom)
            )
            _switch = (Tmp_ss_attributes().o.tmp_ss_size > 6).switch.boolean(
                group_4.o.selection,
                BooleanShrink(Boolean=group_4.o.selection, Shrink=1),
            )
            curve_to_mesh = (
                g.SetHandleType(curve=set_curve_radius)
                >> g.SetSplineResolution(resolution=subdivisions / 4.0)
                >> g.CurveToMesh(
                    profile_curve=g.CurveCircle(
                        resolution=profile_resolution * 2.0, radius=0.31
                    ),
                    scale=named_attribute.o.exists.switch.float(
                        1.0, named_attribute.o.attribute
                    ),
                    fill_caps=True,
                )
            )
        _evaluate_closure = g.EvaluateClosure()
        join_geometry = g.JoinGeometry(geometry=(group_1, curve_to_mesh))
        viewer = g.Viewer()
        group_5 >> viewer

        join_geometry >> geometry


class MN_utils_style_cartoon(CustomGeometryGroup):
    _name = ".MN_utils_style_cartoon"
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry._mn_utils_style_cartoon"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        shade_smooth = tree.inputs.boolean(
            "Shade Smooth",
            True,
            description="Apply smooth shading to the created geometry",
        )
        _interpolate_color = tree.inputs.boolean(
            "Interpolate Color",
            True,
            description="Interpolate between distinct color selections",
        )
        material = tree.inputs.material(
            "Material", description="Material to apply to the resulting geometry"
        )
        ca_curve_threshold = tree.inputs.float(
            "CA Curve Threshold",
            4.5,
            description="Distance (Angstroms) over which subsequent CA points are treated as a new chain",
            min_value=0.0,
            max_value=10_000.0,
        )
        with tree.inputs.panel("Profile"):
            _profile_curve = tree.inputs.geometry(
                "Profile Curve", description="A custom curve-cirlce making SS ribbons."
            )
            profile_resolution = tree.inputs.integer(
                "Profile Resolution", 4, min_value=4, max_value=100
            )
        with tree.inputs.panel("Cylinder"):
            as_cylinders = tree.inputs.boolean("As Cylinders", False)
            _cylinder_curved = tree.inputs.boolean("Cylinder Curved", True)
            _cylinder_radius = tree.inputs.float(
                "Cylinder Radius", 2.0, min_value=0.0, max_value=10_000.0
            )
            _cylinder_resolution = tree.inputs.integer(
                "Cylinder Resolution", 12, min_value=3, max_value=512
            )
            _cylinder_subdivisions = tree.inputs.integer(
                "Cylinder Subdivisions", 5, min_value=1
            )
        with tree.inputs.panel("Helix"):
            helix_thickness = tree.inputs.float(
                "Helix Thickness", 0.5, min_value=0.0, max_value=10_000.0
            )
            helix_width = tree.inputs.float("Helix Width", 2.0)
            helix_subdivisions = tree.inputs.integer(
                "Helix Subdivisions", 5, min_value=1, max_value=20
            )
            _helix_smoothing = tree.inputs.boolean(
                "Helix smoothing",
                True,
                description="Smoothen out AH to be more cylindrical.",
            )
        with tree.inputs.panel("Arrows"):
            as_arrows = tree.inputs.boolean(
                "As Arrows",
                False,
                description="Render beta-strands with directional arrows.",
            )
            arrows_sharp = tree.inputs.menu("Arrows Sharp")
            _arrows_point = tree.inputs.boolean("Arrows Point", False)
            _arrow_thickness_scale = tree.inputs.float(
                "Arrow Thickness Scale", 1.0, min_value=0.0, max_value=10_000.0
            )
            _arrow_width_scale = tree.inputs.float(
                "Arrow Width Scale", 1.0, min_value=-10_000.0, max_value=10_000.0
            )
        with tree.inputs.panel("Sheet"):
            _sheet_rotate = tree.inputs.float("Sheet Rotate", 0.0)
            sheet_thickness = tree.inputs.float("Sheet Thickness", 0.5, min_value=0.0)
            sheet_width = tree.inputs.float(
                "Sheet Width", 2.0, min_value=0.0, max_value=10_000.0
            )
            sheet_smoothing = tree.inputs.float(
                "Sheet Smoothing", 1.0, min_value=0.0, max_value=1.0
            )
            _sheet_subdivision = tree.inputs.integer(
                "Sheet Subdivision", 3, min_value=1, max_value=20
            )
        with tree.inputs.panel("Loop"):
            loop_radius = tree.inputs.float(
                "Loop Radius", 0.3, min_value=0.0, max_value=3.0
            )
            loop_subdivisions = tree.inputs.integer("Loop Subdivisions", 6, min_value=1)
            loop_resolution = tree.inputs.integer(
                "Loop Resolution", 8, min_value=3, max_value=512
            )
        cartoon_mesh = tree.outputs.geometry("Cartoon Mesh")
        ca_splines = tree.outputs.geometry("CA Splines")

        group = AtomsToCACurves(
            atoms=atoms,
            selection=selection,
            bs_smoothing=sheet_smoothing,
            threshold=ca_curve_threshold,
        )
        group_1 = CAToLoops(
            Geometry=group,
            Subdivisions=loop_subdivisions,
            Radius=loop_radius,
            **{"As Cylinders": as_cylinders, "Profile Resolution": loop_resolution},
        )
        group_2 = SplitCurves(Curves=group, Selection=IsHelix().o.selection, Expand=0)
        group_3 = NodeGroup(
            Curves=SplitCurves(Curves=group, Selection=IsLoop().o.selection, Expand=0),
            **{"Reset Curve Normal": True},
            Rotation=(0.0, 0.0, math.tau / 9),
            Menu="Round",
            Radius=0.05,
        )
        _join_geometry = g.JoinGeometry(
            geometry=(
                NodeGroup(Curves=group_2, Scale=(0.8, 2.5, 1.0)),
                NodeGroup(
                    Curves=SplitCurves(
                        Curves=group, Selection=IsSheet().o.selection, Expand=0
                    )
                ),
                group_3,
            )
        )
        menu_switch = g.MenuSwitch.boolean(
            arrows_sharp, {"Sharp": False, "Round": True}
        )
        group_4 = CAToSheet(
            Curve=group,
            **{
                "Profile Resolution": menu_switch.o.output.switch.integer(
                    4, profile_resolution
                )
            },
            Thickness=sheet_thickness,
            Width=sheet_width,
            Subdivisions=helix_subdivisions,
            Rounded=menu_switch.o.output,
            Arrows=as_arrows,
        )
        group_5 = CAToHelix(
            Curve=group,
            Boolean=as_cylinders,
            Thickness=helix_thickness,
            Width=helix_width,
            **{
                "Profile Resolution": (
                    as_cylinders | menu_switch.o.output
                ).switch.integer(4, profile_resolution)
            },
            Subdivisions=helix_subdivisions,
        )
        store_named_attribute = g.StoreNamedAttribute.edge.boolean(
            g.JoinGeometry(geometry=(group_5, group_4, group_1)),
            name="sharp_edge",
            value=~shade_smooth | (g.EdgeAngle().o.signed_angle > math.pi / 3),
        )
        (
            Cleanup(
                geometry=store_named_attribute,
                color_source=group,
                material=material,
                shade_smooth=shade_smooth,
            )
            >> cartoon_mesh
        )
        viewer = g.Viewer()
        group_2 >> viewer

        group >> ca_splines

        arrows_sharp.default_value = "Sharp"


class StyleCartoon(AssetGeometryGroup):
    """
    Style Cartoon

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this style to, discarding unselected points
    quality : InputInteger
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    peptide_shape : InputMenu | Literal["Sharp", "Round"]
        Create rounded sheets and helices
    helix_shape : InputMenu | Literal["Spiral", "Cylinder"]
        Use cylinders for helices instead of ribbons
    helix_thickness : InputFloat
        Thickness for the sheets and helices
    helix_width : InputFloat
        Width for the sheets and helices
    sheet_arrows : InputMenu | Literal["Arrow", "Flat"]
        User arrows for sheets
    sheet_thickness : InputFloat
        Sheet Thickness
    sheet_width : InputFloat
        Sheet Width
    sheet_smoothing : InputFloat
        Smoothing to apply to sheets
    loop_radius : InputFloat
        Radius of the loops for unstructure regions
    backbone_shape : InputMenu | Literal["Cylinder", "Rectangle"]
        Backbone Shape
    nucleic_width : InputFloat
        Nucleic Width
    nucleic_thickness : InputFloat
        Nucleic Thickness
    nucleic_radius : InputFloat
        Nucleic Radius
    base_shape : InputMenu | Literal["Cylinder", "Rectangle"]
        Base Shape
    base_scale_cylinder : InputVector
        Base Scale Cylinder
    base_scale_rectangle : InputVector
        Base Scale Rectangle
    color_blur : InputBoolean
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : InputBoolean
        Apply smooth shading to the created geometry
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this style to, discarding unselected points
    i.quality : IntegerSocket
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    i.peptide_shape : MenuSocket
        Create rounded sheets and helices
    i.helix_shape : MenuSocket
        Use cylinders for helices instead of ribbons
    i.helix_thickness : FloatSocket
        Thickness for the sheets and helices
    i.helix_width : FloatSocket
        Width for the sheets and helices
    i.sheet_arrows : MenuSocket
        User arrows for sheets
    i.sheet_thickness : FloatSocket
        Sheet Thickness
    i.sheet_width : FloatSocket
        Sheet Width
    i.sheet_smoothing : FloatSocket
        Smoothing to apply to sheets
    i.loop_radius : FloatSocket
        Radius of the loops for unstructure regions
    i.backbone_shape : MenuSocket
        Backbone Shape
    i.nucleic_width : FloatSocket
        Nucleic Width
    i.nucleic_thickness : FloatSocket
        Nucleic Thickness
    i.nucleic_radius : FloatSocket
        Nucleic Radius
    i.base_shape : MenuSocket
        Base Shape
    i.base_scale_cylinder : VectorSocket
        Base Scale Cylinder
    i.base_scale_rectangle : VectorSocket
        Base Scale Rectangle
    i.color_blur : BooleanSocket
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        The generated geometry for the style node group
    """

    _name = "Style Cartoon"
    _asset_name = "Style Cartoon"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "node_tool_idname": "geometry.style_cartoon",
        "show_modifier_manage_panel": False,
        "is_modifier": True,
    }

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this style to, discarding unselected points"""
        quality: IntegerSocket
        """A lower value results in less geometry, with a higher value meaning better looking but more dense geometry"""
        peptide_shape: MenuSocket
        """Create rounded sheets and helices"""
        helix_shape: MenuSocket
        """Use cylinders for helices instead of ribbons"""
        helix_thickness: FloatSocket
        """Thickness for the sheets and helices"""
        helix_width: FloatSocket
        """Width for the sheets and helices"""
        sheet_arrows: MenuSocket
        """User arrows for sheets"""
        sheet_thickness: FloatSocket
        """Sheet Thickness"""
        sheet_width: FloatSocket
        """Sheet Width"""
        sheet_smoothing: FloatSocket
        """Smoothing to apply to sheets"""
        loop_radius: FloatSocket
        """Radius of the loops for unstructure regions"""
        backbone_shape: MenuSocket
        """Backbone Shape"""
        nucleic_width: FloatSocket
        """Nucleic Width"""
        nucleic_thickness: FloatSocket
        """Nucleic Thickness"""
        nucleic_radius: FloatSocket
        """Nucleic Radius"""
        base_shape: MenuSocket
        """Base Shape"""
        base_scale_cylinder: VectorSocket
        """Base Scale Cylinder"""
        base_scale_rectangle: VectorSocket
        """Base Scale Rectangle"""
        color_blur: BooleanSocket
        """Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating"""
        shade_smooth: BooleanSocket
        """Apply smooth shading to the created geometry"""
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
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        quality: InputInteger = 2,
        peptide_shape: InputMenu | Literal["Sharp", "Round"] = "Sharp",
        helix_shape: InputMenu | Literal["Spiral", "Cylinder"] = "Spiral",
        helix_thickness: InputFloat = 0.6,
        helix_width: InputFloat = 2.2,
        sheet_arrows: InputMenu | Literal["Arrow", "Flat"] = "Arrow",
        sheet_thickness: InputFloat = 0.5,
        sheet_width: InputFloat = 2.5,
        sheet_smoothing: InputFloat = 0.5,
        loop_radius: InputFloat = 0.4,
        backbone_shape: InputMenu | Literal["Cylinder", "Rectangle"] = "Cylinder",
        nucleic_width: InputFloat = 3.0,
        nucleic_thickness: InputFloat = 1.0,
        nucleic_radius: InputFloat = 1.5,
        base_shape: InputMenu | Literal["Cylinder", "Rectangle"] = "Rectangle",
        base_scale_cylinder: InputVector = None,
        base_scale_rectangle: InputVector = None,
        color_blur: InputBoolean = False,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Quality": quality,
                "Peptide Shape": peptide_shape,
                "Helix Shape": helix_shape,
                "Helix Thickness": helix_thickness,
                "Helix Width": helix_width,
                "Sheet Arrows": sheet_arrows,
                "Sheet Thickness": sheet_thickness,
                "Sheet Width": sheet_width,
                "Sheet Smoothing": sheet_smoothing,
                "Loop Radius": loop_radius,
                "Backbone Shape": backbone_shape,
                "Nucleic Width": nucleic_width,
                "Nucleic Thickness": nucleic_thickness,
                "Nucleic Radius": nucleic_radius,
                "Base Shape": base_shape,
                "Base Scale Cylinder": base_scale_cylinder,
                "Base Scale Rectangle": base_scale_rectangle,
                "Color Blur": color_blur,
                "Shade Smooth": shade_smooth,
                "Material": material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this style to, discarding unselected points",
            hide_value=True,
        )
        quality = tree.inputs.integer(
            "Quality",
            2,
            description="A lower value results in less geometry, with a higher value meaning better looking but more dense geometry",
            min_value=0,
            max_value=8,
        )
        with tree.inputs.panel("Peptide", default_closed=True):
            peptide_shape = tree.inputs.menu(
                "Peptide Shape",
                description="Create rounded sheets and helices",
                expanded=True,
                optional_label=True,
            )
            with tree.inputs.panel("Helix", default_closed=True):
                helix_shape = tree.inputs.menu(
                    "Helix Shape",
                    description="Use cylinders for helices instead of ribbons",
                    expanded=True,
                    optional_label=True,
                )
                helix_thickness = tree.inputs.float(
                    "Helix Thickness",
                    0.6,
                    description="Thickness for the sheets and helices",
                    min_value=0.0,
                )
                helix_width = tree.inputs.float(
                    "Helix Width",
                    2.2,
                    description="Width for the sheets and helices",
                    min_value=0.0,
                )
            with tree.inputs.panel("Sheet", default_closed=True):
                sheet_arrows = tree.inputs.menu(
                    "Sheet Arrows",
                    description="User arrows for sheets",
                    expanded=True,
                    optional_label=True,
                )
                sheet_thickness = tree.inputs.float(
                    "Sheet Thickness", 0.5, min_value=0.0
                )
                sheet_width = tree.inputs.float(
                    "Sheet Width", 2.5, min_value=0.0, max_value=10_000.0
                )
                sheet_smoothing = tree.inputs.float(
                    "Sheet Smoothing",
                    0.5,
                    description="Smoothing to apply to sheets",
                    min_value=0.0,
                    max_value=1.0,
                    subtype="FACTOR",
                )
            with tree.inputs.panel("Loop", default_closed=True):
                loop_radius = tree.inputs.float(
                    "Loop Radius",
                    0.4,
                    description="Radius of the loops for unstructure regions",
                    min_value=0.0,
                    max_value=3.0,
                )
        with tree.inputs.panel("Nucleic", default_closed=True):
            backbone_shape = tree.inputs.menu(
                "Backbone Shape", expanded=True, optional_label=True
            )
            nucleic_width = tree.inputs.float(
                "Nucleic Width", 3.0, min_value=0.0, max_value=10_000.0
            )
            nucleic_thickness = tree.inputs.float(
                "Nucleic Thickness", 1.0, min_value=0.0, max_value=10_000.0
            )
            nucleic_radius = tree.inputs.float(
                "Nucleic Radius", 1.5, min_value=0.0, subtype="DISTANCE"
            )
            with tree.inputs.panel("Base", default_closed=True):
                base_shape = tree.inputs.menu(
                    "Base Shape", expanded=True, optional_label=True
                )
                base_scale_cylinder = tree.inputs.vector(
                    "Base Scale Cylinder", (1.0, 1.0, 7.0)
                )
                base_scale_rectangle = tree.inputs.vector(
                    "Base Scale Rectangle", (2.5, 0.8, 7.0)
                )
        with tree.inputs.panel("Material", default_closed=True):
            color_blur = tree.inputs.boolean(
                "Color Blur",
                False,
                description="Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating",
            )
            shade_smooth = tree.inputs.boolean(
                "Shade Smooth",
                True,
                description="Apply smooth shading to the created geometry",
                structure_type="SINGLE",
                force_non_field=True,
            )
            material = tree.inputs.material(
                "Material",
                description="Material to apply to the resulting geometry",
                optional_label=True,
            )
        geometry = tree.outputs.geometry(
            "Geometry", description="The generated geometry for the style node group"
        )

        _capture = g.CaptureAttribute.point()
        closure_zone = g.ClosureZone()
        atoms_1 = closure_zone.inputs.geometry("Atoms")
        geometry_1 = closure_zone.outputs.geometry("Geometry")
        capture_1 = g.CaptureAttribute.point(geometry=atoms_1, selection=selection)
        group = SeparatePolymers(atoms=capture_1.o.geometry)
        math_1 = quality * 3.0
        math_2 = quality * 5.0
        menu_switch = g.MenuSwitch.integer(base_shape, {"Cylinder": 0, "Rectangle": 1})
        menu_switch_1 = g.MenuSwitch.integer(
            backbone_shape, {"Cylinder": 0, "Rectangle": 1}
        )
        index_switch = g.IndexSwitch.vector(
            menu_switch_1.o.output,
            (
                (0.0, 0.0, 0.0),
                g.CombineXYZ(y=nucleic_width, z=nucleic_thickness, x=1.0),
            ),
        )
        group_1 = MN_utils_style_ribbon_nucleic(
            atoms=group.o.nucleic,
            selection=capture_1.o.selection,
            material=material,
            switch=g.IndexSwitch.boolean(menu_switch_1.o.output, (True, False)),
            backbone_subdivisions=quality * 2.0,
            backbone_resolution=quality * 4.0,
            backbone_radius=g.IndexSwitch.float(
                menu_switch_1.o.output, (nucleic_radius, 0.0)
            ),
            backbone_shade_smooth=shade_smooth,
            backbone_scale=index_switch,
            base_scale=g.IndexSwitch.vector(
                menu_switch.o.output, (base_scale_cylinder, base_scale_rectangle)
            ),
            base_resolution=g.IndexSwitch.integer(menu_switch.o.output, (12, 4)),
        )
        group_2 = MN_utils_style_cartoon(
            Atoms=group.o.peptide,
            Selection=capture_1.o.selection,
            **{"Shade Smooth": shade_smooth, "Interpolate Color": color_blur},
            Material=material,
            **{
                "Profile Resolution": quality * 4.0,
                "As Cylinders": g.MenuSwitch.boolean(
                    helix_shape, {"Spiral": False, "Cylinder": True}
                ).o.output,
                "Cylinder Radius": helix_width,
                "Cylinder Resolution": math_2,
                "Cylinder Subdivisions": math_1,
                "Helix Thickness": helix_thickness,
                "Helix Width": helix_width,
                "Helix Subdivisions": math_1,
                "As Arrows": g.MenuSwitch.boolean(
                    sheet_arrows, {"Arrow": True, "Flat": False}
                ).o.output,
                "Arrows Sharp": peptide_shape,
                "Arrow Thickness Scale": 1.31,
                "Arrow Width Scale": 1.03,
                "Sheet Thickness": sheet_thickness,
                "Sheet Width": sheet_width,
                "Sheet Smoothing": sheet_smoothing,
                "Sheet Subdivision": math_1,
                "Loop Radius": loop_radius,
                "Loop Subdivisions": math_1,
                "Loop Resolution": math_2,
            },
        )
        (
            g.JoinGeometry(geometry=(group_2.o.cartoon_mesh, group_1.o.geometry))
            >> geometry_1
        )
        EvaluateOnAtoms(geometry=atoms, closure=closure_zone.closure) >> geometry

        peptide_shape.default_value = "Sharp"
        helix_shape.default_value = "Spiral"
        sheet_arrows.default_value = "Arrow"
        backbone_shape.default_value = "Cylinder"
        base_shape.default_value = "Rectangle"


ASSET = StyleCartoon

ASSET_METADATA = {
    "catalog_id": "541e6649-2ea6-4225-b1ee-5c0da6f5f1f6",
}
