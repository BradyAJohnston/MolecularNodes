# Node-group asset 'DNA From Curve' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING, Literal
import bpy
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputGeometry, InputInteger, InputMenu
from .angstrom_to_world import AngstromToWorld
from .color_backbone import ColorBackbone
from .color_res_name import ColorResName
from .curve_rotation import CurveRotation
from .curve_visualize import CurveVisualize
from .edge_length import EdgeLength
from .is_hydrogen import IsHydrogen
from .is_side_chain import IsSideChain
from .offset_curve import OffsetCurve
from .residue_id import ResidueID
from .set_color import SetColor
from .transform_local_axis import TransformLocalAxis


class DNASequenceToID(CustomGeometryGroup):
    _name = "DNA Sequence to ID"

    def _build_group(self, tree):
        string = tree.inputs.string("String", "", optional_label=True)
        length = tree.outputs.integer("Length")
        res_id = tree.outputs.integer(
            "res_id", description="Output list with evaluated field values"
        )

        special_characters = g.SpecialCharacters()
        replace_string = (
            string.uppercase()
            .replace(g.String(string=" "))
            .replace(special_characters.o.line_break)
        )
        trim_string = g.TrimString(
            string=g.TrimString(string=replace_string.replace(special_characters.o.tab))
        )
        replace_string_1 = (
            trim_string.o.string.replace("A").replace("C").replace("G").replace("T")
        )
        format_string = g.String(
            string="Non-standard base letters detected: {s}"
        ).o.string.format({"s": replace_string_1})
        slice_string = trim_string.o.string.slice(g.Index(), 1)
        trim_string.o.string.length() >> length
        switch = g.Compare.string.equal(slice_string, "C").o.result.switch.integer(
            g.Compare.string.equal(slice_string, "A").o.result.switch.integer(-1), 1
        )
        (
            g.Compare.string.equal(slice_string, "T").o.result.switch.integer(
                g.Compare.string.equal(slice_string, "G").o.result.switch.integer(
                    switch, 2
                ),
                3,
            )
            >> res_id
        )
        _warning = g.Warning.warning(replace_string_1.length() > 0, format_string)


class CustomWorldObjectSpace(CustomGeometryGroup):
    _name = "Custom/World/Object Space"
    _color_tag = "CONVERTER"
    _tree_properties = {"default_group_node_width": 200}

    def _build_group(self, tree):
        space = tree.inputs.menu("Space", optional_label=True)
        custom_0_world_1_object_2 = tree.outputs.integer(
            "Custom = 0,World = 1, Object = 2"
        )
        is_custom_space = tree.outputs.boolean(
            "Is Custom Space",
            description="True if this item is chosen by the menu input",
        )
        is_world_space = tree.outputs.boolean(
            "Is World Space",
            description="True if this item is chosen by the menu input",
        )
        is_object_space = tree.outputs.boolean(
            "Is Object Space",
            description="True if this item is chosen by the menu input",
        )

        menu_switch = g.MenuSwitch.integer(
            space,
            {
                "Custom Space": (0, "Use a user-defined space."),
                "World Space": (1, "Use world-space coordinates."),
                "Object Space": (
                    2,
                    "Use the coordinates in object space of the provided object.",
                ),
            },
        )

        menu_switch >> custom_0_world_1_object_2
        menu_switch.o.custom_space >> is_custom_space
        menu_switch.o.world_space >> is_world_space
        menu_switch.o.object_space >> is_object_space

        space.default_value = "Custom Space"


class CustomForce(CustomGeometryGroup):
    _name = "Custom Force"
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Create a custom force field to be used in simulations."
    }

    def _build_group(self, tree):
        mode = tree.inputs.menu(
            "Mode",
            description="How the force field is defined.",
            optional_label=True,
            structure_type="SINGLE",
            force_non_field=True,
        )
        closure = tree.inputs.closure(
            "Closure",
            description="The closure computing the force at each point of the input geometry.",
            structure_type="SINGLE",
            force_non_field=True,
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of points this force is applied to.",
            hide_value=True,
            structure_type="FIELD",
        )
        force = tree.inputs.vector(
            "Force",
            (0.0, 0.0, 0.0),
            description="The actual force vector.",
            structure_type="FIELD",
            subtype="XYZ",
        )
        with tree.inputs.panel(
            "Geometry Space",
            description="Geometry space settings.",
            default_closed=True,
        ):
            geometry_space = tree.inputs.menu(
                "Geometry Space",
                description="The space the geometry is transformed to before the field is evaluated.",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
            geometry_space_1 = tree.inputs.matrix(
                "Geometry Space",
                description="Custom geometry space.",
                structure_type="SINGLE",
                force_non_field=True,
            )
            object = tree.inputs.object(
                "Object",
                description="Object to take the geometry space from.",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
                default_input="SELF_OBJECT",
            )
        with tree.inputs.panel(
            "Force Space", description="Force space settings.", default_closed=True
        ):
            force_space = tree.inputs.menu(
                "Force Space",
                description="The space the provided force vector is in.",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
            force_space_1 = tree.inputs.matrix(
                "Force Space",
                description="Custom force space.",
                structure_type="SINGLE",
                force_non_field=True,
            )
            object_1 = tree.inputs.object(
                "Object",
                description="Object to take the force space from.",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
                default_input="SELF_OBJECT",
            )
        with tree.inputs.panel(
            "Filter", description="Filter settings.", default_closed=True
        ):
            filter = tree.inputs.string(
                "Filter",
                "",
                description="Comma-separated list of tags this effector should be applied to.",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
        force_1 = tree.outputs.bundle(
            "Force",
            description="A bundle containing all information about this effector.",
        )

        closure_zone = g.ClosureZone()
        geometry = closure_zone.inputs.geometry("Geometry")
        to_world_transform = closure_zone.inputs.matrix("To World Transform")
        geometry_1 = closure_zone.outputs.geometry("Geometry")
        selection_1 = closure_zone.outputs.boolean("Selection")
        force_2 = closure_zone.outputs.vector("Force")
        with g.Frame("Used for type inferencing"):
            evaluate_closure = g.EvaluateClosure(closure, define_signature=True)
            evaluate_closure.inputs.geometry("Geometry", structure_type="SINGLE")
            evaluate_closure.inputs.matrix(
                "To World Transform", structure_type="SINGLE"
            )
            evaluate_closure.outputs.geometry("Geometry", structure_type="SINGLE")
            evaluate_closure.outputs.boolean("Selection", structure_type="FIELD")
            evaluate_closure.outputs.vector("Force", structure_type="FIELD")
        with g.Frame("World space to geometry space"):
            index_switch = g.IndexSwitch.matrix(
                CustomWorldObjectSpace(
                    Space=geometry_space
                ).o.custom_0_world_1_object_2,
                (
                    geometry_space_1.invert(),
                    None,
                    g.ObjectInfo(object=object).o.transform,
                ),
            )
            invert_matrix = index_switch.o.output.invert()
            _string = g.String(
                string="The geometry nodes to be transformed into the space requested by the user before evaluating the fields."
            )
        with g.Frame("Force space to world space"):
            index_switch_1 = g.IndexSwitch.matrix(
                CustomWorldObjectSpace(Space=force_space).o.custom_0_world_1_object_2,
                (force_space_1, None, g.ObjectInfo(object=object_1).o.transform),
            )
            _string_1 = g.String(
                string="The computed force field needs to be transformed back into the space of the simulation."
            )
        transform_direction = g.MultiplyMatrices(
            matrix=to_world_transform.invert(), matrix_001=index_switch_1
        ).o.matrix.transform_direction(force)
        transform_geometry = geometry >> g.TransformGeometry(
            transform=g.MultiplyMatrices(
                matrix=invert_matrix, matrix_001=to_world_transform
            ),
            mode="Matrix",
        )
        transform_geometry >> geometry_1
        selection >> selection_1
        transform_direction >> force_2
        menu_switch = g.MenuSwitch.closure(
            mode,
            {
                "Field": (
                    closure_zone.closure,
                    "Create the force just by providing a field and specifying the relevant spaces.",
                ),
                "Closure": (
                    closure,
                    "Create a force by providing a closure that computes it.",
                ),
            },
        )
        combine_bundle = g.CombineBundle()
        combine_bundle.items.string("Type", "Blender.Force")
        combine_bundle.items.string("filter", filter)
        combine_bundle.items.closure("closure", menu_switch.o.output)

        combine_bundle.o.bundle >> force_1

        mode.default_value = "Field"
        geometry_space.default_value = "World Space"
        force_space.default_value = "World Space"


class CurveSegment(CustomGeometryGroup):
    _name = "Curve Segment"
    _color_tag = "INPUT"

    def _build_group(self, tree):
        segment_length = tree.outputs.float(
            "Segment Length", description="Distance to previous point on curve"
        )
        segment_direction = tree.outputs.vector(
            "Segment Direction",
            description="Direction from previous neighboring point on segment",
        )
        neighbor_index = tree.outputs.integer(
            "Neighbor Index",
            description="Index of previous neighboring point on segment",
        )

        evaluate_on_domain = g.Boolean(boolean=True).o.boolean.spline.evaluate()
        offset_point_in_curve = g.OffsetPointInCurve(offset=-1)
        position = g.Position()
        switch = offset_point_in_curve.o.is_valid_offset.switch.integer(
            g.Index(), offset_point_in_curve.o.point_index
        )
        (
            evaluate_on_domain.switch.integer(true=switch.point.evaluate())
            >> neighbor_index
        )
        evaluate_on_domain_1 = (
            position.o.position - position.o.position.point.at(switch)
        ).point.evaluate()
        (
            evaluate_on_domain.switch.float(true=evaluate_on_domain_1.length())
            >> segment_length
        )
        (
            evaluate_on_domain.switch.vector(
                (0.0, 0.0, 0.0), evaluate_on_domain_1.normalize()
            )
            >> segment_direction
        )


class StoreSegmentLength(CustomGeometryGroup):
    _name = "Store Segment Length"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree):
        curves = tree.inputs.geometry("Curves")
        name = tree.inputs.string("Name", "rest_length", optional_label=True)
        curves_1 = tree.outputs.geometry("Curves")

        (
            curves
            >> g.StoreNamedAttribute.point.float(
                name=name, value=CurveSegment().o.segment_length
            )
            >> curves_1
        )


class StoreEdgeLength(CustomGeometryGroup):
    _name = "Store Edge Length"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree):
        curves = tree.inputs.geometry("Curves")
        name = tree.inputs.string("Name", "rest_length", optional_label=True)
        curves_1 = tree.outputs.geometry("Curves")

        (
            curves
            >> g.StoreNamedAttribute.edge.float(name=name, value=EdgeLength())
            >> curves_1
        )


class StoreSegmentRotation(CustomGeometryGroup):
    _name = "Store Segment Rotation"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree):
        curves = tree.inputs.geometry("Curves")
        name = tree.inputs.string("Name", "rest_rotation", optional_label=True)
        curves_1 = tree.outputs.geometry("Curves")

        position = g.Position()
        axes_to_rotation = g.AxesToRotation(
            primary_axis=position.o.position.point.at(
                g.OffsetPointInCurve(offset=1).o.point_index
            )
            - position,
            secondary_axis=g.Normal().o.normal,
        )
        (
            curves
            >> g.StoreNamedAttribute.point.quaternion(name=name, value=axes_to_rotation)
            >> curves_1
        )


class StoreBendRotation(CustomGeometryGroup):
    _name = "Store Bend Rotation"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree):
        curves = tree.inputs.geometry("Curves")
        bend_rotation_name = tree.inputs.string(
            "Bend Rotation Name", "rest_bend_rotation", optional_label=True
        )
        rotation_name = tree.inputs.string(
            "Rotation Name", "rest_rotation", optional_label=True
        )
        curves_1 = tree.outputs.geometry("Curves")

        named_attribute = g.NamedAttribute.quaternion(rotation_name)
        rotate_rotation = named_attribute.o.attribute.invert().rotate(
            named_attribute.o.attribute.point.at(
                g.OffsetPointInCurve(offset=1).o.point_index
            ),
            rotation_space="LOCAL",
        )
        (
            curves
            >> g.StoreNamedAttribute.point.quaternion(
                name=bend_rotation_name, value=rotate_rotation
            )
            >> curves_1
        )


class SetupStructuralRestData(CustomGeometryGroup):
    _name = "Setup Structural Rest Data"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 180}

    def _build_group(self, tree):
        geometry = tree.inputs.geometry("Geometry")
        geometry_1 = tree.outputs.geometry("Geometry")

        get_geometry_component = g.GetGeometryComponent(geometry=geometry, type="Curve")
        get_geometry_component_1 = get_geometry_component >> g.GetGeometryComponent(
            type="Grease Pencil"
        )
        get_geometry_component_2 = get_geometry_component_1 >> g.GetGeometryComponent()
        group = StoreSegmentRotation(
            Curves=StoreSegmentLength(Curves=get_geometry_component.o.component)
        )
        group_1 = StoreSegmentRotation(
            Curves=StoreSegmentLength(Curves=get_geometry_component_1.o.component)
        )
        join_geometry = g.JoinGeometry(
            geometry=(
                get_geometry_component_2.o.geometry,
                StoreEdgeLength(Curves=get_geometry_component_2.o.component),
                StoreBendRotation(Curves=group_1),
                StoreBendRotation(Curves=group),
            )
        )

        join_geometry >> geometry_1


class Damping(CustomGeometryGroup):
    _name = "Damping"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree):
        linear = tree.inputs.float("Linear", 2.0, min_value=0.0, structure_type="FIELD")
        angular = tree.inputs.float(
            "Angular", 2.0, min_value=0.0, structure_type="FIELD"
        )
        with tree.inputs.panel("Filter", default_closed=True):
            filter = tree.inputs.string(
                "Filter", "", structure_type="SINGLE", force_non_field=True
            )
        damping = tree.outputs.bundle("Damping")

        combine_bundle = g.CombineBundle()
        combine_bundle.items.string("Type", "Blender.Damping")
        combine_bundle.items.string("filter", filter)
        combine_bundle.items.float("linear_damping", linear)
        combine_bundle.items.float("angular_damping", angular)

        combine_bundle.o.bundle >> damping


class RodBendTwistConstraint(CustomGeometryGroup):
    _name = "Rod Bend/Twist Constraint"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 200}

    def _build_group(self, tree):
        compliance = tree.inputs.float(
            "Compliance", 0.0001, min_value=0.0, structure_type="FIELD"
        )
        with tree.inputs.panel("Custom Bend Rotation", default_closed=True):
            custom_bend_rotation = tree.inputs.boolean(
                "Custom Bend Rotation",
                False,
                structure_type="FIELD",
                is_panel_toggle=True,
            )
            bend_rotation = tree.inputs.rotation(
                "Bend Rotation", (0.0, 0.0, 0.0), structure_type="FIELD"
            )
        with tree.inputs.panel("Filter", default_closed=True):
            filter = tree.inputs.string(
                "Filter",
                "",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
        with tree.inputs.panel("Debug", default_closed=True):
            error_threshold = tree.inputs.float(
                "Error Threshold",
                0.01,
                min_value=0.00001,
                structure_type="SINGLE",
                subtype="ANGLE",
                force_non_field=True,
            )
        constraint = tree.outputs.bundle("Constraint")

        with g.Frame("Rest rotation"):
            switch = custom_bend_rotation.switch.rotation(
                g.NamedAttribute.quaternion("rest_bend_rotation").o.attribute,
                bend_rotation,
            )
        combine_bundle = g.CombineBundle()
        combine_bundle.items.string("Type", "Blender.Constraint.RodBendTwist")
        combine_bundle.items.string("filter", filter)
        combine_bundle.items.rotation("rest_bend_rotation", switch)
        combine_bundle.items.float("compliance", compliance)
        combine_bundle.items.float("error_threshold", error_threshold)

        combine_bundle.o.bundle >> constraint


class RodStretchShearConstraint(CustomGeometryGroup):
    _name = "Rod Stretch/Shear Constraint"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 200}

    def _build_group(self, tree):
        compliance = tree.inputs.float(
            "Compliance", 0.0001, min_value=0.0, structure_type="FIELD"
        )
        with tree.inputs.panel("Custom Length", default_closed=True):
            custom_length = tree.inputs.boolean(
                "Custom Length", False, is_panel_toggle=True
            )
            length = tree.inputs.float("Length", 0.0, structure_type="FIELD")
        with tree.inputs.panel("Filter", default_closed=True):
            filter = tree.inputs.string(
                "Filter",
                "",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
        with tree.inputs.panel("Debug", default_closed=True):
            error_threshold = tree.inputs.float(
                "Error Threshold",
                0.001,
                min_value=0.00001,
                structure_type="SINGLE",
                subtype="DISTANCE",
                force_non_field=True,
            )
            position_lambda = tree.inputs.string(
                "Position Lambda",
                "",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
            rotation_lambda = tree.inputs.string(
                "Rotation Lambda",
                "",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
        constraint = tree.outputs.bundle("Constraint")

        with g.Frame("Rest length"):
            switch = custom_length.switch.float(
                g.NamedAttribute.float("rest_length").o.attribute, length
            )
        combine_bundle = g.CombineBundle()
        combine_bundle.items.string("Type", "Blender.Constraint.RodStretchShear")
        combine_bundle.items.string("filter", filter)
        combine_bundle.items.float("rest_length", switch)
        combine_bundle.items.float("compliance", compliance)
        combine_bundle.items.float("error_threshold", error_threshold)
        combine_bundle.items.string("lambda_position_attribute", position_lambda)
        combine_bundle.items.string("lambda_rotation_attribute", rotation_lambda)

        combine_bundle.o.bundle >> constraint


class SoftnessToCompliance(CustomGeometryGroup):
    _name = "Softness to Compliance"
    _color_tag = "CONVERTER"
    _tree_properties = {"default_group_node_width": 160}

    def _build_group(self, tree):
        softness = tree.inputs.float(
            "Softness", 0.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        compliance = tree.outputs.float("Compliance")

        softness.max(0.0) ** 5.0 >> compliance


class SetEffector(CustomGeometryGroup):
    _name = "Set Effector"
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Add effector information to the geometry bundle.",
        "default_group_node_width": 160,
    }

    def _build_group(self, tree):
        geometry = tree.inputs.geometry(
            "Geometry",
            description="Geometry to add the effector to.",
            structure_type="SINGLE",
            force_non_field=True,
        )
        effector = tree.inputs.bundle(
            "Effector",
            description="The actual bundle containing one more multiple effectors.",
            structure_type="SINGLE",
            force_non_field=True,
        )
        name = tree.inputs.string(
            "Name",
            "",
            description="Name of the effector in the geometry bundle.",
            optional_label=True,
            structure_type="SINGLE",
            force_non_field=True,
        )
        geometry_1 = tree.outputs.geometry(
            "Geometry",
            description="The same geometry but with additional effector information.",
        )

        get_geometry_bundle = g.GetGeometryBundle(geometry=geometry, remove=True)
        join_strings = g.JoinStrings(
            (
                g.String(string="effectors"),
                g.Compare.string.equal(name, "").o.result.switch.string(
                    name, "default"
                ),
            ),
            delimiter="/",
        )
        (
            get_geometry_bundle
            >> g.SetGeometryBundle(
                bundle=g.StoreBundleItem.bundle(
                    get_geometry_bundle.o.bundle, join_strings, effector
                )
            )
            >> geometry_1
        )


class Collider(CustomGeometryGroup):
    _name = "Collider"
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Turn a geometry into a collider.",
        "default_group_node_width": 200,
        "show_modifier_manage_panel": False,
        "is_modifier": True,
    }

    def _build_group(self, tree):
        geometry = tree.inputs.geometry(
            "Geometry",
            description="The collider geometry.",
            structure_type="SINGLE",
            force_non_field=True,
        )
        deforming = tree.inputs.boolean(
            "Deforming",
            False,
            description="Whether the collider geometry is changing.",
            structure_type="SINGLE",
            force_non_field=True,
        )
        boundary = tree.inputs.boolean(
            "Boundary",
            False,
            description="If true, the simulated geometry should stay within the mesh instead of it being pushed out of it.",
            structure_type="SINGLE",
            force_non_field=True,
        )
        edge_contacts = tree.inputs.boolean(
            "Edge Contacts",
            False,
            description="Increased collision accuracy (only for collisions with curves/hair currently).",
            structure_type="SINGLE",
            force_non_field=True,
        )
        margin = tree.inputs.float(
            "Margin",
            0.0,
            description="Additional margin around the collider to prevent penetration.",
            min_value=0.0,
            structure_type="SINGLE",
            subtype="DISTANCE",
            force_non_field=True,
        )
        friction = tree.inputs.float(
            "Friction",
            0.2,
            description="Friction coefficient for this collider.",
            min_value=0.0,
            max_value=1.0,
            structure_type="SINGLE",
            subtype="FACTOR",
            force_non_field=True,
        )
        softness = tree.inputs.float(
            "Softness",
            0.0,
            description="Increasing softness allows some penetration.",
            min_value=0.0,
            max_value=1.0,
            structure_type="SINGLE",
            subtype="FACTOR",
            force_non_field=True,
        )
        with tree.inputs.panel(
            "Filter", description="Filter settings.", default_closed=True
        ):
            filter = tree.inputs.string(
                "Filter",
                "",
                description="Comma-separated list of tags this collider should apply to.",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
        with tree.inputs.panel(
            "Debug", description="Debug settings.", default_closed=True
        ):
            error_threshold = tree.inputs.float(
                "Error Threshold",
                0.001,
                description="Used to compute a relative error of this constraint to determine solve quality.",
                min_value=0.00001,
                structure_type="SINGLE",
                subtype="DISTANCE",
                force_non_field=True,
            )
        geometry_1 = tree.outputs.geometry(
            "Geometry",
            description="The input geometry with additional collider information attached.",
        )
        collider = tree.outputs.bundle(
            "Collider",
            description="The collider effector bundle that can be passed into a simulation.",
        )

        with g.Frame("Prepare collider geometry"):
            _string = g.String(
                string="By convention, collider geometries are always in world space."
            )
            transform_geometry = g.GeometryToInstance(geometry) >> g.TransformGeometry(
                transform=g.ObjectInfo(object=g.SelfObject()).o.transform, mode="Matrix"
            )
        combine_bundle = g.CombineBundle()
        combine_bundle.items.string("Type", "Blender.Collider.Mesh")
        combine_bundle.items.string("filter", filter)
        combine_bundle.items.geometry("geometry", transform_geometry)
        combine_bundle.items.float("margin", margin)
        combine_bundle.items.float("friction", friction)
        combine_bundle.items.float(
            "compliance", SoftnessToCompliance(Softness=softness)
        )
        combine_bundle.items.boolean("deforming", deforming)
        combine_bundle.items.boolean("use_edge_contacts", edge_contacts)
        combine_bundle.items.boolean("is_boundary", boundary)
        combine_bundle.items.float("error_threshold", error_threshold)
        (
            SetEffector(
                Geometry=geometry, Effector=combine_bundle.o.bundle, Name="collider"
            )
            >> geometry_1
        )

        combine_bundle.o.bundle >> collider


class ObjectEffector(CustomGeometryGroup):
    _name = "Object Effector"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree):
        object = tree.inputs.object("Object", optional_label=True)
        effector = tree.outputs.bundle("Effector")

        get_bundle_item = g.GetBundleItem.bundle(
            g.GetGeometryBundle(
                geometry=g.ObjectInfo(object=object).o.geometry
            ).o.bundle,
            "effectors",
        )

        get_bundle_item.o.item >> effector


class CollectionEffector(CustomGeometryGroup):
    _name = "Collection Effector"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree):
        collection = tree.inputs.collection("Collection", optional_label=True)
        with tree.inputs.panel("Naming", default_closed=True):
            name_pattern = tree.inputs.string(
                "Name Pattern", "effector_{}", optional_label=True
            )
        effectors = tree.outputs.bundle("Effectors")

        collection_children = g.CollectionChildren(
            collection=collection, recursive=True
        )
        repeat_zone = g.RepeatZone(collection_children.o.objects.list_length())
        bundle = repeat_zone.items.bundle("Bundle")
        store_bundle_item = g.StoreBundleItem.bundle(
            bundle.current,
            name_pattern.format({"i": repeat_zone.iteration}),
            ObjectEffector(Object=collection_children.o.objects[repeat_zone.iteration]),
        )
        store_bundle_item >> bundle.next

        bundle.result >> effectors


class PinPositions(CustomGeometryGroup):
    _name = "Pin Positions"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 200}

    def _build_group(self, tree):
        selection = tree.inputs.boolean(
            "Selection", True, hide_value=True, structure_type="FIELD"
        )
        position = tree.inputs.vector(
            "Position", (0.0, 0.0, 0.0), structure_type="FIELD"
        )
        compliance = tree.inputs.float(
            "Compliance", 0.0, min_value=0.0, structure_type="FIELD"
        )
        with tree.inputs.panel("Filter", default_closed=True):
            filter = tree.inputs.string(
                "Filter",
                "",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
        with tree.inputs.panel("Debug", default_closed=True):
            error_threshold = tree.inputs.float(
                "Error Threshold",
                0.001,
                min_value=0.00001,
                structure_type="SINGLE",
                subtype="DISTANCE",
                force_non_field=True,
            )
            lambda_ = tree.inputs.string(
                "Lambda",
                "",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
        constraint = tree.outputs.bundle("Constraint")

        combine_bundle = g.CombineBundle()
        combine_bundle.items.string("Type", "Blender.Constraint.PinPosition")
        combine_bundle.items.string("filter", filter)
        combine_bundle.items.boolean("selection", selection)
        combine_bundle.items.vector("position", position)
        combine_bundle.items.float("compliance", compliance)
        combine_bundle.items.float("error_threshold", error_threshold)
        combine_bundle.items.string("lambda_attribute", lambda_)

        combine_bundle.o.bundle >> constraint


class PinRotation(CustomGeometryGroup):
    _name = "Pin Rotation"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 200}

    def _build_group(self, tree):
        selection = tree.inputs.boolean("Selection", True, structure_type="FIELD")
        rotation = tree.inputs.rotation(
            "Rotation", (0.0, 0.0, 0.0), structure_type="FIELD"
        )
        compliance = tree.inputs.float(
            "Compliance", 0.0, min_value=0.0, structure_type="FIELD"
        )
        with tree.inputs.panel("Filter", default_closed=True):
            filter = tree.inputs.string(
                "Filter",
                "",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
        with tree.inputs.panel("Debug", default_closed=True):
            error_threshold = tree.inputs.float(
                "Error Threshold", 0.01, min_value=0.00001, subtype="ANGLE"
            )
        constraint = tree.outputs.bundle("Constraint")

        combine_bundle = g.CombineBundle()
        combine_bundle.items.string("Type", "Blender.Constraint.PinRotation")
        combine_bundle.items.string("filter", filter)
        combine_bundle.items.boolean("selection", selection)
        combine_bundle.items.rotation("rotation", rotation)
        combine_bundle.items.float("compliance", compliance)
        combine_bundle.items.float("error_threshold", error_threshold)

        combine_bundle.o.bundle >> constraint


class StringToList(CustomGeometryGroup):
    _name = "String to List"
    _color_tag = "CONVERTER"

    def _build_group(self, tree):
        string = tree.inputs.string("String", "", optional_label=True)
        separator = tree.inputs.string("Separator", ",", optional_label=True)
        list = tree.outputs.string("List")

        trim_string = g.TrimString(
            string=g.SplitString(string=string, separator=separator)
        )

        trim_string >> list


class SetGeometryTags(CustomGeometryGroup):
    _name = "Set Geometry Tags"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree):
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to get the bundle of"
        )
        tags = tree.inputs.string("Tags", "", optional_label=True)
        geometry_1 = tree.outputs.geometry("Geometry")

        get_geometry_bundle = g.GetGeometryBundle(geometry=geometry, remove=True)
        store_bundle_item = g.StoreBundleItem.string(
            get_geometry_bundle.o.bundle, "tags", StringToList(String=tags)
        )
        (
            get_geometry_bundle
            >> g.SetGeometryBundle(bundle=store_bundle_item)
            >> geometry_1
        )


class SetFriction(CustomGeometryGroup):
    _name = "Set Friction"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 160}

    def _build_group(self, tree):
        geometry = tree.inputs.geometry(
            "Geometry",
            description="Geometry to store a new attribute with the given name on",
        )
        static = tree.inputs.float("Static", 0.5, min_value=0.0)
        dynamic = tree.inputs.float("Dynamic", 0.5, min_value=0.0)
        geometry_1 = tree.outputs.geometry("Geometry")

        (
            geometry
            >> g.StoreNamedAttribute.point.float(name="static_friction", value=static)
            >> g.StoreNamedAttribute.point.float(name="dynamic_friction", value=dynamic)
            >> geometry_1
        )


class ThinRodMomentOfInertia(CustomGeometryGroup):
    _name = "Thin Rod Moment of Inertia"
    _color_tag = "CONVERTER"

    def _build_group(self, tree):
        mass = tree.inputs.float("Mass", 0.0, hide_value=True)
        moment_of_inertia = tree.outputs.vector("Moment of Inertia")

        offset_point_in_curve = g.OffsetPointInCurve(offset=1)
        position = g.Position()
        vector_math = (
            position.o.position.point.at(offset_point_in_curve.o.point_index) - position
        )
        math_1 = (
            mass + mass.point.at(offset_point_in_curve.o.point_index)
        ) * vector_math.dot(vector_math)
        math_2 = math_1 * 0.08333334
        combine_xyz = g.CombineXYZ(x=math_2, y=math_2)

        combine_xyz >> moment_of_inertia


class SetMass(CustomGeometryGroup):
    _name = "Set Mass"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 160}

    def _build_group(self, tree):
        geometry = tree.inputs.geometry("Geometry")
        mass = tree.inputs.float(
            "Mass", 1.0, min_value=0.0, structure_type="FIELD", subtype="MASS"
        )
        with tree.inputs.panel("Moment of Inertia", default_closed=True):
            moment_of_inertia = tree.inputs.boolean(
                "Moment of Inertia",
                False,
                structure_type="SINGLE",
                is_panel_toggle=True,
                force_non_field=True,
            )
            moment_of_inertia_mode = tree.inputs.menu(
                "Moment of Inertia Mode", optional_label=True
            )
            moment_of_inertia_1 = tree.inputs.vector(
                "Moment of Inertia", (1.0, 1.0, 1.0), optional_label=True
            )
        geometry_1 = tree.outputs.geometry("Geometry")

        with g.Frame("Simple Rod Model for Moments of Inertia"):
            _named_attribute = g.NamedAttribute.float("mass")
        menu_switch = g.MenuSwitch.vector(
            moment_of_inertia_mode,
            {
                "Custom": moment_of_inertia_1,
                "Thin Rod": ThinRodMomentOfInertia(
                    Mass=g.NamedAttribute.float("mass").o.attribute
                ),
            },
        )
        store_named_attribute = g.StoreNamedAttribute.point.float(
            geometry, name="mass", value=mass
        )
        store_named_attribute_1 = g.StoreNamedAttribute.point.vector(
            store_named_attribute, name="moment_of_inertia", value=menu_switch.o.output
        )
        (
            moment_of_inertia.switch.geometry(
                store_named_attribute, store_named_attribute_1
            )
            >> geometry_1
        )

        moment_of_inertia_mode.default_value = "Custom"


class SimAttributes(CustomGeometryGroup):
    _name = "Sim Attributes"
    _color_tag = "CONVERTER"

    def _build_group(self, tree):
        extra = tree.inputs.string("Extra", "")
        sim_attributes = tree.outputs.string("Sim Attributes")

        join_strings = g.JoinStrings(
            (
                g.String(
                    string="position, rotation, velocity, angular_velocity, sim:*"
                ),
                extra,
            ),
            delimiter=",",
        )
        StringToList(String=join_strings) >> sim_attributes


class ApplyGeoCacheDeformOnly(CustomGeometryGroup):
    _name = "Apply Geo Cache - Deform Only"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 240}

    def _build_group(self, tree):
        new_geometry = tree.inputs.geometry("New Geometry")
        cache_geometry = tree.inputs.geometry("Cache Geometry")
        sim_attributes = tree.inputs.string(
            "Sim Attributes",
            "",
            description="List of attribute names (not) to transfer. A wildcard (*) at the end is allowed",
            optional_label=True,
        )
        post_sim = tree.inputs.boolean("Post Sim", False)
        mixed_geometry = tree.outputs.geometry("Mixed Geometry")

        transfer_attributes = g.TransferAttributes(
            target=new_geometry, source=cache_geometry, attribute_names=sim_attributes
        )
        transfer_attributes.node.warning_propagation = "ERRORS_AND_WARNINGS"
        transfer_attributes_1 = g.TransferAttributes(
            target=new_geometry,
            source=cache_geometry,
            attribute_names=g.FieldToList(items={"Value": "*"}),
        )
        transfer_attributes_1.node.warning_propagation = "ERRORS_AND_WARNINGS"
        (
            post_sim.switch.geometry(transfer_attributes, transfer_attributes_1)
            >> mixed_geometry
        )


class SetGeoUpdaterDeformOnly(CustomGeometryGroup):
    _name = "Set Geo Updater - Deform Only"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 220}

    def _build_group(self, tree):
        geometry = tree.inputs.geometry("Geometry")
        extra_sim_attributes = tree.inputs.string(
            "Extra Sim Attributes", "", optional_label=True
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        closure_zone = g.ClosureZone()
        geometry_2 = closure_zone.inputs.geometry("Geometry")
        cache = closure_zone.inputs.geometry("Cache")
        post_sim = closure_zone.inputs.boolean("Post Sim")
        geometry_3 = closure_zone.outputs.geometry("Geometry")
        get_geometry_bundle = g.GetGeometryBundle(geometry=geometry, remove=True)
        group = ApplyGeoCacheDeformOnly(
            **{
                "New Geometry": geometry_2,
                "Cache Geometry": cache,
                "Sim Attributes": SimAttributes(Extra=extra_sim_attributes),
                "Post Sim": post_sim,
            }
        )
        group >> geometry_3
        store_bundle_item = g.StoreBundleItem(
            bundle=get_geometry_bundle.o.bundle,
            item=closure_zone.closure,
            path="sim_apply_cache",
            socket_type="CLOSURE",
            structure_type="SINGLE",
        )
        (
            get_geometry_bundle
            >> g.SetGeometryBundle(bundle=store_bundle_item)
            >> geometry_1
        )


class IsEffectorForGeometry(CustomGeometryGroup):
    _name = "Is Effector for Geometry"
    _color_tag = "CONVERTER"
    _tree_properties = {"default_group_node_width": 180}

    def _build_group(self, tree):
        effector = tree.inputs.bundle("Effector")
        _effector_path = tree.inputs.string(
            "Effector Path",
            "",
            optional_label=True,
            structure_type="SINGLE",
            force_non_field=True,
        )
        geometry = tree.inputs.geometry("Geometry")
        _geometry_path = tree.inputs.string(
            "Geometry Path",
            "",
            optional_label=True,
            structure_type="SINGLE",
            force_non_field=True,
        )
        affects_geometry = tree.outputs.boolean("Affects Geometry")

        get_bundle_item = g.GetBundleItem.string(
            g.GetGeometryBundle(geometry=geometry).o.bundle, "tags"
        )
        get_bundle_item.node.warning_propagation = "NONE"
        separate_bundle = g.SeparateBundle(effector)
        filter = separate_bundle.items.string("filter")
        tag_filter = g.TagFilter(tag_filter=filter, tags=get_bundle_item.o.item)

        tag_filter >> affects_geometry


class ForEachSimGeometry(CustomGeometryGroup):
    _name = "For Each Sim Geometry"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 160}

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        closure = tree.inputs.closure("Closure")
        world_1 = tree.outputs.bundle("World")

        get_nested_bundle_paths = g.GetNestedBundlePaths(
            bundle=world, mode="Data Type", data_type="Geometry"
        )
        repeat_zone = g.RepeatZone(get_nested_bundle_paths.o.paths.list_length())
        world_2 = repeat_zone.items.bundle("World", world)
        get_list_item = get_nested_bundle_paths.o.paths[repeat_zone.iteration]
        get_bundle_item = g.GetBundleItem.geometry(world_2.current, get_list_item, True)
        evaluate_closure = g.EvaluateClosure(closure, define_signature=True)
        evaluate_closure.inputs.bundle(
            "World", get_bundle_item.o.bundle, structure_type="SINGLE"
        )
        evaluate_closure.inputs.geometry(
            "Geometry", get_bundle_item.o.item, structure_type="SINGLE"
        )
        evaluate_closure.inputs.string("Path", get_list_item, structure_type="SINGLE")
        world_3 = evaluate_closure.outputs.bundle("World", structure_type="SINGLE")
        geometry = evaluate_closure.outputs.geometry(
            "Geometry", structure_type="SINGLE"
        )
        g.StoreBundleItem.geometry(world_3, get_list_item, geometry) >> world_2.next

        world_2.result >> world_1


class EvaluateCustomEffectors(CustomGeometryGroup):
    _name = "Evaluate Custom Effectors"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree):
        world = tree.inputs.bundle(
            "World", structure_type="SINGLE", force_non_field=True
        )
        stage = tree.inputs.string(
            "Stage",
            "",
            optional_label=True,
            structure_type="SINGLE",
            force_non_field=True,
        )
        to_world_transform = tree.inputs.matrix(
            "To World Transform", structure_type="SINGLE", force_non_field=True
        )
        world_1 = tree.outputs.bundle("World")

        with g.Frame("Evaluate Geometry Effectors"):
            closure_zone = g.ClosureZone()
            world_2 = closure_zone.inputs.bundle("World", structure_type="SINGLE")
            geometry = closure_zone.inputs.geometry("Geometry", structure_type="SINGLE")
            path = closure_zone.inputs.string("Path", structure_type="SINGLE")
            world_3 = closure_zone.outputs.bundle("World", structure_type="SINGLE")
            geometry_1 = closure_zone.outputs.geometry(
                "Geometry", structure_type="SINGLE"
            )
            get_nested_bundle_paths = g.GetNestedBundlePaths(
                bundle=world,
                mode="Bundle Type",
                bundle_type="Blender.CustomEffector.Geometry",
            )
            repeat_zone = g.RepeatZone(get_nested_bundle_paths.o.paths.list_length())
            world_4 = repeat_zone.items.bundle("World", world)
            get_list_item = g.GetListItem(
                list=get_nested_bundle_paths,
                index=repeat_zone.iteration,
                socket_type="STRING",
                structure_type="SINGLE",
            )
            get_bundle_item = g.GetBundleItem.bundle(world_4.current, get_list_item)
            group = IsEffectorForGeometry(
                Effector=get_bundle_item.o.item,
                **{"Effector Path": get_list_item},
                Geometry=geometry,
                **{"Geometry Path": path},
            )
            separate_bundle = g.SeparateBundle(get_bundle_item.o.item)
            closure = separate_bundle.items.closure("closure")
            stage_1 = separate_bundle.items.string("stage")
            evaluate_closure = g.EvaluateClosure(closure)
            evaluate_closure.inputs.geometry(
                "Geometry", geometry, structure_type="SINGLE"
            )
            evaluate_closure.inputs.matrix(
                "To World Transform", to_world_transform, structure_type="SINGLE"
            )
            geometry_2 = evaluate_closure.outputs.geometry(
                "Geometry", structure_type="SINGLE"
            )
            world_2 >> world_3
            group.o.affects_geometry.switch.geometry(geometry, geometry_2) >> geometry_1
            switch = g.Compare.string.equal(stage_1, stage).o.result.switch.bundle(
                world_4.current,
                ForEachSimGeometry(World=world_4.current, Closure=closure_zone.closure),
            )
            switch >> world_4.next
        with g.Frame("Evaluate World Effectors"):
            get_nested_bundle_paths_1 = g.GetNestedBundlePaths(
                bundle=world_4.result,
                mode="Bundle Type",
                bundle_type="Blender.CustomEffector.World",
            )
            repeat_zone_1 = g.RepeatZone(
                get_nested_bundle_paths_1.o.paths.list_length()
            )
            world_5 = repeat_zone_1.items.bundle("World", world_4.result)
            get_list_item_1 = g.GetListItem(
                list=get_nested_bundle_paths_1,
                index=repeat_zone_1.iteration,
                socket_type="STRING",
                structure_type="SINGLE",
            )
            separate_bundle_1 = g.SeparateBundle(
                g.GetBundleItem.bundle(world_5.current, get_list_item_1).o.item
            )
            stage_2 = separate_bundle_1.items.string("stage", structure_type="SINGLE")
            closure_1 = separate_bundle_1.items.closure("closure")
            evaluate_closure_1 = g.EvaluateClosure(closure_1)
            evaluate_closure_1.inputs.bundle(
                "World", world_5.current, structure_type="SINGLE"
            )
            evaluate_closure_1.inputs.matrix("To World Transform", to_world_transform)
            world_6 = evaluate_closure_1.outputs.bundle(
                "World", structure_type="SINGLE"
            )
            (
                g.Compare.string.equal(stage, stage_2).o.result.switch.bundle(
                    world_5.current, world_6
                )
                >> world_5.next
            )

        world_5.result >> world_1


class RenameSimAttributes(CustomGeometryGroup):
    _name = "Rename Sim Attributes"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 160}

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        mode = tree.inputs.menu("Mode", optional_label=True)
        old = tree.inputs.string("Old", "", optional_label=True)
        new = tree.inputs.string("New", "", optional_label=True)
        world_1 = tree.outputs.bundle("World")

        closure_zone = g.ClosureZone()
        world_2 = closure_zone.inputs.bundle("World", structure_type="SINGLE")
        geometry = closure_zone.inputs.geometry("Geometry", structure_type="SINGLE")
        closure_zone.inputs.string("Path", structure_type="SINGLE")
        world_3 = closure_zone.outputs.bundle("World", structure_type="SINGLE")
        geometry_1 = closure_zone.outputs.geometry("Geometry", structure_type="SINGLE")
        world_2 >> world_3
        (
            g.RenameAttribute(
                geometry=geometry, mode=mode, old=old, new=new, overwrite=True
            )
            >> geometry_1
        )
        ForEachSimGeometry(World=world, Closure=closure_zone.closure) >> world_1

        mode.default_value = "Single"


class SetPreviousWorldItems(CustomGeometryGroup):
    _name = "Set Previous World Items"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 180}

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        cache = tree.inputs.bundle("Cache")
        world_1 = tree.outputs.bundle("World")

        get_nested_bundle_paths = g.GetNestedBundlePaths(
            bundle=world,
            mode="Bundle Type",
            pattern_mode="Wildcard",
            bundle_type="*",
            data_type="Bundle",
        )
        repeat_zone = g.RepeatZone(get_nested_bundle_paths.o.paths.list_length())
        world_2 = repeat_zone.items.bundle("World", world)
        get_list_item = get_nested_bundle_paths.o.paths[repeat_zone.iteration]
        get_bundle_item = g.GetBundleItem.bundle(world_2.current, get_list_item, True)
        get_bundle_item_1 = g.GetBundleItem.bundle(cache, get_list_item)
        store_bundle_item = g.StoreBundleItem.bundle(
            world_2.current,
            get_list_item,
            g.StoreBundleItem.bundle(
                get_bundle_item.o.item, "previous", get_bundle_item_1.o.item
            ),
        )
        switch = (get_bundle_item.o.exists & get_bundle_item_1.o.exists).switch.bundle(
            world_2.current, store_bundle_item
        )
        switch >> world_2.next

        world_2.result >> world_1


class ApplyEachGeometryCache(CustomGeometryGroup):
    _name = "Apply Each Geometry Cache"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 200}

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        cache = tree.inputs.bundle("Cache")
        post_sim = tree.inputs.boolean("Post Sim", False)
        world_1 = tree.outputs.bundle("World")

        closure_zone = g.ClosureZone()
        world_2 = closure_zone.inputs.bundle("World", structure_type="SINGLE")
        geometry = closure_zone.inputs.geometry("Geometry", structure_type="SINGLE")
        path = closure_zone.inputs.string("Path", structure_type="SINGLE")
        world_3 = closure_zone.outputs.bundle("World", structure_type="SINGLE")
        geometry_1 = closure_zone.outputs.geometry("Geometry", structure_type="SINGLE")
        get_geometry_bundle = g.GetGeometryBundle(geometry=geometry)
        get_bundle_item = g.GetBundleItem.geometry(cache, path)
        evaluate_closure = g.EvaluateClosure(
            g.GetBundleItem.closure(
                get_geometry_bundle.o.bundle, "sim_apply_cache"
            ).o.item
        )
        evaluate_closure.inputs.geometry("Geometry", get_geometry_bundle.o.geometry)
        evaluate_closure.inputs.geometry("Cache", get_bundle_item.o.item)
        evaluate_closure.inputs.boolean("Post Sim", post_sim)
        geometry_2 = evaluate_closure.outputs.geometry("Geometry")
        world_2 >> world_3
        (
            get_bundle_item.o.exists.switch.geometry(get_geometry_bundle, geometry_2)
            >> geometry_1
        )
        ForEachSimGeometry(World=world, Closure=closure_zone.closure) >> world_1


class EvaluateEffectorAttributes(CustomGeometryGroup):
    _name = "Evaluate Effector Attributes"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 180}

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        world_1 = tree.outputs.bundle("World")

        with g.Frame("Pin Position Constraint"):
            closure_zone = g.ClosureZone()
            geometry = closure_zone.inputs.geometry("Geometry")
            effector = closure_zone.inputs.bundle("Effector")
            effector_path = closure_zone.inputs.string("Effector Path")
            geometry_1 = closure_zone.outputs.geometry("Geometry")
            separate_bundle = g.SeparateBundle(effector)
            selection = separate_bundle.items.boolean("selection")
            position = separate_bundle.items.vector("position")
            compliance = separate_bundle.items.float("compliance")
            capture = g.CaptureAttribute.point(geometry=geometry, selection=selection)
            position_1 = capture.items.vector("position", position)
            compliance_1 = capture.items.float("compliance", compliance)
            store_named_attribute = (
                capture.o.geometry
                >> g.StoreNamedAttribute.point.boolean(
                    name=g.FormatString(
                        "sim:prop:{}:selection", items={"p": effector_path}
                    ),
                    value=capture.o.selection,
                )
                >> g.StoreNamedAttribute.point.vector(
                    name=g.FormatString(
                        "sim:prop:{}:position", items={"p": effector_path}
                    ),
                    value=position_1.output,
                )
                >> g.StoreNamedAttribute.point.float(
                    name=g.FormatString(
                        "sim:prop:{}:compliance", items={"p": effector_path}
                    ),
                    value=compliance_1.output,
                )
            )
            store_named_attribute >> geometry_1
            combine_bundle = g.CombineBundle()
            combine_bundle.items.string("Type", "Blender.FieldPreEvaluation")
            combine_bundle.items.string(
                "effector_type", g.String(string="Blender.Constraint.PinPosition")
            )
            combine_bundle.items.closure("closure", closure_zone.closure)
        with g.Frame("Pin Rotation Constraint"):
            closure_zone_1 = g.ClosureZone()
            geometry_2 = closure_zone_1.inputs.geometry("Geometry")
            effector_1 = closure_zone_1.inputs.bundle("Effector")
            effector_path_1 = closure_zone_1.inputs.string("Effector Path")
            geometry_3 = closure_zone_1.outputs.geometry("Geometry")
            separate_bundle_1 = g.SeparateBundle(effector_1)
            selection_1 = separate_bundle_1.items.boolean("selection")
            rotation = separate_bundle_1.items.rotation("rotation")
            compliance_2 = separate_bundle_1.items.float("compliance")
            capture_1 = g.CaptureAttribute.point(
                geometry=geometry_2, selection=selection_1
            )
            position_2 = capture_1.items.rotation("position", rotation)
            compliance_3 = capture_1.items.float("compliance", compliance_2)
            store_named_attribute_1 = (
                capture_1.o.geometry
                >> g.StoreNamedAttribute.point.boolean(
                    name=g.FormatString(
                        "sim:prop:{}:selection", items={"p": effector_path_1}
                    ),
                    value=capture_1.o.selection,
                )
                >> g.StoreNamedAttribute.point.quaternion(
                    name=g.FormatString(
                        "sim:prop:{}:rotation", items={"p": effector_path_1}
                    ),
                    value=position_2.output,
                )
                >> g.StoreNamedAttribute.point.float(
                    name=g.FormatString(
                        "sim:prop:{}:compliance", items={"p": effector_path_1}
                    ),
                    value=compliance_3.output,
                )
            )
            store_named_attribute_1 >> geometry_3
            combine_bundle_1 = g.CombineBundle()
            combine_bundle_1.items.string("Type", "Blender.FieldPreEvaluation")
            combine_bundle_1.items.string(
                "effector_type", g.String(string="Blender.Constraint.PinRotation")
            )
            combine_bundle_1.items.closure("closure", closure_zone_1.closure)
        with g.Frame("Damping"):
            closure_zone_2 = g.ClosureZone()
            geometry_4 = closure_zone_2.inputs.geometry("Geometry")
            effector_2 = closure_zone_2.inputs.bundle("Effector")
            effector_path_2 = closure_zone_2.inputs.string("Effector Path")
            geometry_5 = closure_zone_2.outputs.geometry("Geometry")
            separate_bundle_2 = g.SeparateBundle(effector_2)
            linear_damping = separate_bundle_2.items.float("linear_damping")
            angular_damping = separate_bundle_2.items.float("angular_damping")
            capture_2 = g.CaptureAttribute.point(geometry=geometry_4)
            linear_damping_1 = capture_2.items.float("linear_damping", linear_damping)
            angular_damping_1 = capture_2.items.float(
                "angular_damping", angular_damping
            )
            store_named_attribute_2 = (
                capture_2.o.geometry
                >> g.StoreNamedAttribute.point.float(
                    name=g.FormatString(
                        "sim:prop:{}:linear", items={"p": effector_path_2}
                    ),
                    value=linear_damping_1.output,
                )
                >> g.StoreNamedAttribute.point.float(
                    name=g.FormatString(
                        "sim:prop:{}:angular", items={"p": effector_path_2}
                    ),
                    value=angular_damping_1.output,
                )
            )
            store_named_attribute_2 >> geometry_5
            combine_bundle_2 = g.CombineBundle()
            combine_bundle_2.items.string("Type", "Blender.FieldPreEvaluation")
            combine_bundle_2.items.string(
                "effector_type", g.String(string="Blender.Damping")
            )
            combine_bundle_2.items.closure("closure", closure_zone_2.closure)
        with g.Frame("Edge Length Constraint"):
            closure_zone_3 = g.ClosureZone()
            geometry_6 = closure_zone_3.inputs.geometry("Geometry")
            effector_3 = closure_zone_3.inputs.bundle("Effector")
            effector_path_3 = closure_zone_3.inputs.string("Effector Path")
            geometry_7 = closure_zone_3.outputs.geometry("Geometry")
            separate_bundle_3 = g.SeparateBundle(effector_3)
            rest_length = separate_bundle_3.items.float("rest_length")
            compliance_4 = separate_bundle_3.items.float("compliance")
            capture_3 = g.CaptureAttribute.edge(geometry=geometry_6)
            rest_length_1 = capture_3.items.float("rest_length", rest_length)
            compliance_5 = capture_3.items.float("compliance", compliance_4)
            store_named_attribute_3 = (
                capture_3.o.geometry
                >> g.StoreNamedAttribute.edge.float(
                    name=g.FormatString(
                        "sim:prop:{}:rest_length", items={"p": effector_path_3}
                    ),
                    value=rest_length_1.output,
                )
                >> g.StoreNamedAttribute.edge.float(
                    name=g.FormatString(
                        "sim:prop:{}:compliance", items={"p": effector_path_3}
                    ),
                    value=compliance_5.output,
                )
            )
            store_named_attribute_3 >> geometry_7
            combine_bundle_3 = g.CombineBundle()
            combine_bundle_3.items.string("Type", "Blender.FieldPreEvaluation")
            combine_bundle_3.items.string(
                "effector_type", g.String(string="Blender.Constraint.EdgeLength")
            )
            combine_bundle_3.items.closure("closure", closure_zone_3.closure)
        with g.Frame("Rod Stretch Shear Constraint"):
            closure_zone_4 = g.ClosureZone()
            geometry_8 = closure_zone_4.inputs.geometry("Geometry")
            effector_4 = closure_zone_4.inputs.bundle("Effector")
            effector_path_4 = closure_zone_4.inputs.string("Effector Path")
            geometry_9 = closure_zone_4.outputs.geometry("Geometry")
            separate_bundle_4 = g.SeparateBundle(effector_4)
            rest_length_2 = separate_bundle_4.items.float("rest_length")
            compliance_6 = separate_bundle_4.items.float("compliance")
            capture_4 = g.CaptureAttribute.point(geometry=geometry_8)
            rest_length_3 = capture_4.items.float("rest_length", rest_length_2)
            compliance_7 = capture_4.items.float("compliance", compliance_6)
            store_named_attribute_4 = (
                capture_4.o.geometry
                >> g.StoreNamedAttribute.point.float(
                    name=g.FormatString(
                        "sim:prop:{}:rest_length", items={"p": effector_path_4}
                    ),
                    value=rest_length_3.output,
                )
                >> g.StoreNamedAttribute.point.float(
                    name=g.FormatString(
                        "sim:prop:{}:compliance", items={"p": effector_path_4}
                    ),
                    value=compliance_7.output,
                )
            )
            store_named_attribute_4 >> geometry_9
            combine_bundle_4 = g.CombineBundle()
            combine_bundle_4.items.string("Type", "Blender.FieldPreEvaluation")
            combine_bundle_4.items.string(
                "effector_type", g.String(string="Blender.Constraint.RodStretchShear")
            )
            combine_bundle_4.items.closure("closure", closure_zone_4.closure)
        with g.Frame("Rod Bend Twist Constraint"):
            closure_zone_5 = g.ClosureZone()
            geometry_10 = closure_zone_5.inputs.geometry("Geometry")
            effector_5 = closure_zone_5.inputs.bundle("Effector")
            effector_path_5 = closure_zone_5.inputs.string("Effector Path")
            geometry_11 = closure_zone_5.outputs.geometry("Geometry")
            separate_bundle_5 = g.SeparateBundle(effector_5)
            rest_bend_rotation = separate_bundle_5.items.rotation("rest_bend_rotation")
            compliance_8 = separate_bundle_5.items.float("compliance")
            capture_5 = g.CaptureAttribute.point(geometry=geometry_10)
            rest_bend_rotation_1 = capture_5.items.rotation(
                "rest_bend_rotation", rest_bend_rotation
            )
            compliance_9 = capture_5.items.float("compliance", compliance_8)
            store_named_attribute_5 = (
                capture_5.o.geometry
                >> g.StoreNamedAttribute.point.quaternion(
                    name=g.FormatString(
                        "sim:prop:{}:rest_bend_rotation", items={"p": effector_path_5}
                    ),
                    value=rest_bend_rotation_1.output,
                )
                >> g.StoreNamedAttribute.point.float(
                    name=g.FormatString(
                        "sim:prop:{}:compliance", items={"p": effector_path_5}
                    ),
                    value=compliance_9.output,
                )
            )
            store_named_attribute_5 >> geometry_11
            combine_bundle_5 = g.CombineBundle()
            combine_bundle_5.items.string("Type", "Blender.FieldPreEvaluation")
            combine_bundle_5.items.string(
                "effector_type", g.String(string="Blender.Constraint.RodBendTwist")
            )
            combine_bundle_5.items.closure("closure", closure_zone_5.closure)
        with g.Frame("Cross Edge Length Constraint"):
            closure_zone_6 = g.ClosureZone()
            geometry_12 = closure_zone_6.inputs.geometry("Geometry")
            effector_6 = closure_zone_6.inputs.bundle("Effector")
            effector_path_6 = closure_zone_6.inputs.string("Effector Path")
            geometry_13 = closure_zone_6.outputs.geometry("Geometry")
            separate_bundle_6 = g.SeparateBundle(effector_6)
            rest_position = separate_bundle_6.items.vector("rest_position")
            compliance_10 = separate_bundle_6.items.float("compliance")
            capture_6 = g.CaptureAttribute.point(geometry=geometry_12)
            rest_position_1 = capture_6.items.vector("rest_position", rest_position)
            compliance_11 = capture_6.items.float("compliance", compliance_10)
            store_named_attribute_6 = (
                capture_6.o.geometry
                >> g.StoreNamedAttribute.point.vector(
                    name=g.FormatString(
                        "sim:prop:{}:rest_position", items={"p": effector_path_6}
                    ),
                    value=rest_position_1.output,
                )
                >> g.StoreNamedAttribute.edge.float(
                    name=g.FormatString(
                        "sim:prop:{}:compliance", items={"p": effector_path_6}
                    ),
                    value=compliance_11.output,
                )
            )
            store_named_attribute_6 >> geometry_13
            combine_bundle_6 = g.CombineBundle()
            combine_bundle_6.items.string("Type", "Blender.FieldPreEvaluation")
            combine_bundle_6.items.string(
                "effector_type", g.String(string="Blender.Constraint.CrossEdgeLength")
            )
            combine_bundle_6.items.closure("closure", closure_zone_6.closure)
        combine_bundle_7 = g.CombineBundle()
        combine_bundle_7.items.bundle("Pin Position", combine_bundle.o.bundle)
        combine_bundle_7.items.bundle("Pin Rotation", combine_bundle_1.o.bundle)
        combine_bundle_7.items.bundle("Damping", combine_bundle_2.o.bundle)
        combine_bundle_7.items.bundle("Edge Length", combine_bundle_3.o.bundle)
        combine_bundle_7.items.bundle("Rod Stretch Shear", combine_bundle_4.o.bundle)
        combine_bundle_7.items.bundle("Rod Bend Twist", combine_bundle_5.o.bundle)
        combine_bundle_7.items.bundle("Cross Edge Length", combine_bundle_6.o.bundle)
        get_nested_bundle_paths = g.GetNestedBundlePaths(
            bundle=combine_bundle_7.o.bundle,
            mode="Bundle Type",
            bundle_type="Blender.FieldPreEvaluation",
        )
        list_length = get_nested_bundle_paths.o.paths.list_length()
        with g.Frame("For each geomery > for each effector type > for each effector"):
            closure_zone_7 = g.ClosureZone()
            world_2 = closure_zone_7.inputs.bundle("World", structure_type="SINGLE")
            geometry_14 = closure_zone_7.inputs.geometry(
                "Geometry", structure_type="SINGLE"
            )
            path = closure_zone_7.inputs.string("Path", structure_type="SINGLE")
            world_3 = closure_zone_7.outputs.bundle("World", structure_type="SINGLE")
            geometry_15 = closure_zone_7.outputs.geometry(
                "Geometry", structure_type="SINGLE"
            )
            repeat_zone = g.RepeatZone(list_length)
            geometry_16 = repeat_zone.items.geometry("Geometry", geometry_14)
            get_bundle_item = g.GetBundleItem.bundle(
                combine_bundle_7.o.bundle,
                get_nested_bundle_paths.o.paths[repeat_zone.iteration],
            )
            separate_bundle_7 = g.SeparateBundle(get_bundle_item.o.item)
            effector_type = separate_bundle_7.items.string("effector_type")
            closure = separate_bundle_7.items.closure("closure")
            get_nested_bundle_paths_1 = g.GetNestedBundlePaths(
                bundle=world_2, bundle_type=effector_type, mode="Bundle Type"
            )
            repeat_zone_1 = g.RepeatZone(
                get_nested_bundle_paths_1.o.paths.list_length()
            )
            geometry_17 = repeat_zone_1.items.geometry("Geometry", geometry_16.current)
            get_list_item = get_nested_bundle_paths_1.o.paths[repeat_zone_1.iteration]
            get_bundle_item_1 = g.GetBundleItem.bundle(world_2, get_list_item)
            evaluate_closure = g.EvaluateClosure(closure)
            evaluate_closure.inputs.geometry("Geometry", geometry_17.current)
            evaluate_closure.inputs.bundle("Effector", get_bundle_item_1.o.item)
            evaluate_closure.inputs.string("Effector Path", get_list_item)
            geometry_18 = evaluate_closure.outputs.geometry("Geometry")
            group = IsEffectorForGeometry(
                Effector=get_bundle_item_1.o.item,
                **{"Effector Path": get_list_item},
                Geometry=geometry_17.current,
                **{"Geometry Path": path},
            )
            (
                group.o.affects_geometry.switch.geometry(
                    geometry_17.current, geometry_18
                )
                >> geometry_17.next
            )
            geometry_17.result >> geometry_16.next
            world_2 >> world_3
            geometry_16.result >> geometry_15
        ForEachSimGeometry(World=world, Closure=closure_zone_7.closure) >> world_1


class EvaluateForces(CustomGeometryGroup):
    _name = "Evaluate Forces"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 160}

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        to_world_transform = tree.inputs.matrix("To World Transform")
        name = tree.inputs.string("Name", "external_force", optional_label=True)
        world_1 = tree.outputs.bundle("World")

        closure_zone = g.ClosureZone()
        world_2 = closure_zone.inputs.bundle("World", structure_type="SINGLE")
        geometry = closure_zone.inputs.geometry("Geometry", structure_type="SINGLE")
        path = closure_zone.inputs.string("Path", structure_type="SINGLE")
        world_3 = closure_zone.outputs.bundle("World", structure_type="SINGLE")
        geometry_1 = closure_zone.outputs.geometry("Geometry", structure_type="SINGLE")
        get_nested_bundle_paths = g.GetNestedBundlePaths(
            bundle=world, mode="Bundle Type", bundle_type="Blender.Force"
        )
        repeat_zone = g.RepeatZone(get_nested_bundle_paths.o.paths.list_length())
        total_force = repeat_zone.items.vector("Total Force")
        get_list_item = get_nested_bundle_paths.o.paths[repeat_zone.iteration]
        get_bundle_item = g.GetBundleItem.bundle(world_2, get_list_item)
        group = IsEffectorForGeometry(
            Effector=get_bundle_item.o.item,
            **{"Effector Path": get_list_item},
            Geometry=geometry,
            **{"Geometry Path": path},
        )
        separate_bundle = g.SeparateBundle(get_bundle_item.o.item)
        closure = separate_bundle.items.closure("closure")
        evaluate_closure = g.EvaluateClosure(closure)
        evaluate_closure.inputs.geometry("Geometry", geometry)
        evaluate_closure.inputs.matrix("To World Transform", to_world_transform)
        geometry_2 = evaluate_closure.outputs.geometry("Geometry")
        selection = evaluate_closure.outputs.boolean("Selection")
        force = evaluate_closure.outputs.vector("Force")
        capture = g.CaptureAttribute.point(geometry=geometry_2, selection=selection)
        force_1 = capture.items.vector("Force", force)
        sample_index = capture.o.geometry >> g.SampleIndex(
            value=force_1.output, index=g.Index(), data_type="FLOAT_VECTOR"
        )
        switch = group.o.affects_geometry.switch.vector(
            total_force.current, total_force.current + sample_index
        )
        switch >> total_force.next
        world_2 >> world_3
        (
            g.StoreNamedAttribute.point.vector(
                geometry, name=name, value=total_force.result
            )
            >> geometry_1
        )
        ForEachSimGeometry(World=world, Closure=closure_zone.closure) >> world_1


class XPBDSolver(CustomGeometryGroup):
    _name = "XPBD Solver"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree):
        world = tree.inputs.bundle(
            "World", description="World state that is updated by the solver"
        )
        delta_time = tree.inputs.float(
            "Delta Time", 0.04, min_value=0.0, subtype="TIME_ABSOLUTE"
        )
        filter = tree.inputs.string(
            "Filter",
            "",
            description="Filters the geometry sets to process based on their tags",
            optional_label=True,
        )
        simulation_to_world = tree.inputs.matrix("Simulation to World")
        with tree.inputs.panel("Solver", default_closed=True):
            substeps = tree.inputs.integer("Substeps", 10, min_value=1)
            constraint_iterations = tree.inputs.integer(
                "Constraint Iterations", 1, min_value=1
            )
            solver_output_path = tree.inputs.string(
                "Solver Output Path",
                "",
                description="Optional output path in the world bundle for solver data",
                optional_label=True,
            )
        world_1 = tree.outputs.bundle("World")

        with g.Frame("Evaluate substeps one at a time"):
            math_1 = 1.0 / substeps
            repeat_zone = g.RepeatZone(substeps)
            world_2 = repeat_zone.items.bundle("World", world)
            math_2 = repeat_zone.iteration * math_1
            xpbd_solver = g.XpbdSolver(
                world=world_2.current,
                delta_time=delta_time / substeps,
                filter=filter,
                simulation_to_world=simulation_to_world,
                constraint_iterations=constraint_iterations,
                solver_path=solver_output_path,
                begin=math_2,
                end=math_2 + math_1,
                substeps=1,
            )
            xpbd_solver >> world_2.next
        with g.Frame("Evaluate all substeps at once"):
            xpbd_solver_1 = g.XpbdSolver(
                world=world,
                delta_time=delta_time,
                filter=filter,
                simulation_to_world=simulation_to_world,
                substeps=substeps,
                constraint_iterations=constraint_iterations,
                solver_path=solver_output_path,
            )
        with g.Frame("The result of approaches should be equal"):
            switch = g.Switch.bundle(false=xpbd_solver_1, true=world_2.result)

        switch >> world_1


class ForEachTypedBundle(CustomGeometryGroup):
    _name = "For Each Typed Bundle"

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        type = tree.inputs.string("Type", "", optional_label=True)
        closure = tree.inputs.closure("Closure")
        world_1 = tree.outputs.bundle("World")

        get_nested_bundle_paths = g.GetNestedBundlePaths(
            bundle=world, bundle_type=type, mode="Bundle Type"
        )
        repeat_zone = g.RepeatZone(get_nested_bundle_paths.o.paths.list_length())
        world_2 = repeat_zone.items.bundle("World", world)
        get_list_item = get_nested_bundle_paths.o.paths[repeat_zone.iteration]
        get_bundle_item = g.GetBundleItem.bundle(world_2.current, get_list_item, True)
        evaluate_closure = g.EvaluateClosure(closure)
        evaluate_closure.inputs.bundle("World", get_bundle_item.o.bundle)
        evaluate_closure.inputs.bundle("Item", get_bundle_item.o.item)
        evaluate_closure.inputs.string("Path", get_list_item)
        world_3 = evaluate_closure.outputs.bundle("World")
        item = evaluate_closure.outputs.bundle("Item")
        g.StoreBundleItem.bundle(world_3, get_list_item, item) >> world_2.next

        world_2.result >> world_1


class SimplifyCachedColliderInfo(CustomGeometryGroup):
    _name = "Simplify Cached Collider Info"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 200}

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        world_1 = tree.outputs.bundle("World")

        with g.Frame("Delete points when the collider doens't deform"):
            closure_zone = g.ClosureZone()
            world_2 = closure_zone.inputs.bundle("World")
            item = closure_zone.inputs.bundle("Item")
            closure_zone.inputs.string("Path")
            world_3 = closure_zone.outputs.bundle("World")
            item_1 = closure_zone.outputs.bundle("Item")
            separate_bundle = g.SeparateBundle(item)
            deforming = separate_bundle.items.boolean("deforming")
            geometry = separate_bundle.items.geometry("geometry")
            remove_named_attribute = deforming.switch.geometry(
                g.DeleteGeometry.point(geometry),
                g.DeleteGeometry.only_edges_faces(geometry),
            ) >> g.RemoveNamedAttribute(pattern_mode="Wildcard", name="*")
            remove_named_attribute.node.warning_propagation = "NONE"
            combine_bundle = g.CombineBundle()
            combine_bundle.items.geometry("geometry", remove_named_attribute)
            world_2 >> world_3
            combine_bundle.o.bundle >> item_1
            (
                ForEachTypedBundle(
                    World=world,
                    Type="Blender.Collider.Mesh",
                    Closure=closure_zone.closure,
                )
                >> world_1
            )


class ClearPreviousWorldItems(CustomGeometryGroup):
    _name = "Clear Previous World Items"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 180}

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        world_1 = tree.outputs.bundle("World")

        get_nested_bundle_paths = g.GetNestedBundlePaths(
            bundle=world,
            mode="Bundle Type",
            pattern_mode="Wildcard",
            bundle_type="*",
            data_type="Bundle",
        )
        repeat_zone = g.RepeatZone(get_nested_bundle_paths.o.paths.list_length())
        world_2 = repeat_zone.items.bundle("World", world)
        join_strings = g.JoinStrings(
            (
                get_nested_bundle_paths.o.paths[repeat_zone.iteration],
                g.String(string="previous"),
            ),
            delimiter="/",
        )
        get_bundle_item = g.GetBundleItem.bundle(world_2.current, join_strings, True)
        get_bundle_item.node.warning_propagation = "ERRORS"
        get_bundle_item >> world_2.next

        world_2.result >> world_1


class RemoveSimAttributes(CustomGeometryGroup):
    _name = "Remove Sim Attributes"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 160}

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        mode = tree.inputs.menu("Mode", optional_label=True)
        name = tree.inputs.string("Name", "", optional_label=True)
        world_1 = tree.outputs.bundle("World")

        closure_zone = g.ClosureZone()
        world_2 = closure_zone.inputs.bundle("World", structure_type="SINGLE")
        geometry = closure_zone.inputs.geometry("Geometry", structure_type="SINGLE")
        closure_zone.inputs.string("Path", structure_type="SINGLE")
        world_3 = closure_zone.outputs.bundle("World", structure_type="SINGLE")
        geometry_1 = closure_zone.outputs.geometry("Geometry", structure_type="SINGLE")
        world_2 >> world_3
        (
            g.RemoveNamedAttribute(geometry=geometry, pattern_mode=mode, name=name)
            >> geometry_1
        )
        ForEachSimGeometry(World=world, Closure=closure_zone.closure) >> world_1

        mode.default_value = "Exact"


class CopySolverData(CustomGeometryGroup):
    _name = "Copy Solver Data"
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Copy solver output data bundle from simulation cache into a world bundle",
        "default_group_node_width": 180,
    }

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        cache = tree.inputs.bundle("Cache")
        world_1 = tree.outputs.bundle("World")

        get_nested_bundle_paths = g.GetNestedBundlePaths(
            bundle=cache, mode="Bundle Type", bundle_type="Blender.XPBDSolverData"
        )
        repeat_zone = g.RepeatZone(get_nested_bundle_paths.o.paths.list_length())
        world_2 = repeat_zone.items.bundle("World", world)
        get_list_item = get_nested_bundle_paths.o.paths[repeat_zone.iteration]
        store_bundle_item = g.StoreBundleItem.bundle(
            world_2.current,
            get_list_item,
            g.GetBundleItem.bundle(cache, get_list_item).o.item,
        )
        store_bundle_item >> world_2.next

        world_2.result >> world_1


class XPBDSimulation(CustomGeometryGroup):
    _name = "XPBD Simulation"
    _color_tag = "GEOMETRY"
    _tree_properties = {"default_group_node_width": 160}

    def _build_group(self, tree):
        world = tree.inputs.bundle("World")
        substeps = tree.inputs.integer("Substeps", 10, min_value=1)
        constraint_steps = tree.inputs.integer("Constraint Steps", 1, min_value=1)
        simulation_to_world = tree.inputs.matrix("Simulation to World")
        time_scale = tree.inputs.float("Time Scale", 1.0, min_value=0.0, max_value=10.0)
        solver_output_path = tree.inputs.string(
            "Solver Output Path",
            "",
            description="Optional output path in the world bundle for solver data",
            optional_label=True,
        )
        world_1 = tree.outputs.bundle("World")

        group = EvaluateCustomEffectors(
            World=world,
            Stage="PRE_SIMULATION",
            **{"To World Transform": simulation_to_world},
        )
        simulation_zone = g.SimulationZone()
        cache = simulation_zone.items.bundle("Cache")
        group_1 = RenameSimAttributes(
            World=cache.current, Mode="Prefix", Old="sim:prop:", New="sim:prop_prev:"
        )
        group_2 = ApplyEachGeometryCache(
            World=SetPreviousWorldItems(World=group, Cache=group_1), Cache=group_1
        )
        group_3 = EvaluateCustomEffectors(
            World=group_2,
            Stage="PRE_SOLVE",
            **{"To World Transform": simulation_to_world},
        )
        group_4 = EvaluateForces(
            World=EvaluateEffectorAttributes(World=group_3),
            **{"To World Transform": simulation_to_world},
        )
        group_5 = XPBDSolver(
            World=group_4,
            **{
                "Delta Time": simulation_zone.delta_time * time_scale,
                "Simulation to World": simulation_to_world,
            },
            Substeps=substeps,
            **{
                "Constraint Iterations": constraint_steps,
                "Solver Output Path": solver_output_path,
            },
        )
        group_5.node.warning_propagation = "ERRORS_AND_WARNINGS"
        group_6 = EvaluateCustomEffectors(
            World=group_5,
            Stage="POST_SOLVE",
            **{"To World Transform": simulation_to_world},
        )
        group_7 = EvaluateCustomEffectors(
            World=group_6,
            Stage="DEFAULT",
            **{"To World Transform": simulation_to_world},
        )
        group_8 = RemoveSimAttributes(
            World=ClearPreviousWorldItems(
                World=SimplifyCachedColliderInfo(World=group_7)
            ),
            Mode="Wildcard",
            Name="sim:prop_prev:*",
        )
        group_8 >> cache.next
        group_9 = RemoveSimAttributes(
            World=cache.result, Mode="Wildcard", Name="sim:prop:*"
        )
        (
            CopySolverData(
                World=ApplyEachGeometryCache(
                    World=world, Cache=group_9, **{"Post Sim": True}
                ),
                Cache=group_9,
            )
            >> world_1
        )


class ConvertSpaceTransform(CustomGeometryGroup):
    _name = "Convert Space Transform"
    _color_tag = "CONVERTER"
    _tree_properties = {"default_group_node_width": 180}

    def _build_group(self, tree):
        from_space = tree.inputs.menu("From Space", optional_label=True)
        from_object = tree.inputs.object("From Object", optional_label=True)
        to_space = tree.inputs.menu("To Space", optional_label=True)
        to_object = tree.inputs.object("To Object", optional_label=True)
        custom_to_world = tree.inputs.matrix("Custom to World")
        transform = tree.outputs.matrix("Transform")
        inverted = tree.outputs.matrix(
            "Inverted",
            description="The inverted matrix or the identity matrix if the input is not invertible",
        )

        with g.Frame('World to "To" space'):
            index_switch = g.IndexSwitch.matrix(
                CustomWorldObjectSpace(Space=to_space).o.custom_0_world_1_object_2,
                (
                    custom_to_world.invert(),
                    None,
                    g.ObjectInfo(object=to_object).o.transform.invert(),
                ),
            )
        with g.Frame('"From" to world space'):
            index_switch_1 = g.IndexSwitch.matrix(
                CustomWorldObjectSpace(Space=from_space).o.custom_0_world_1_object_2,
                (custom_to_world, None, g.ObjectInfo(object=from_object).o.transform),
            )
        multiply_matrices = g.MultiplyMatrices(
            matrix=index_switch, matrix_001=index_switch_1
        )
        multiply_matrices.o.matrix.invert() >> inverted

        multiply_matrices >> transform

        from_space.default_value = "World Space"
        to_space.default_value = "World Space"


class SimulateDNAGuide(CustomGeometryGroup):
    _name = "Simulate DNA Guide"
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Simulate hair attached to a surface while taking effectors such as forces and colliders into account.",
        "default_group_node_width": 200,
        "is_modifier": True,
    }

    def _build_group(self, tree):
        dna_curve = tree.inputs.geometry(
            "DNA Curve",
            description="Input hair curves.",
            structure_type="SINGLE",
            force_non_field=True,
        )
        pin_rotations = tree.inputs.boolean("Pin Rotations", False)
        pin_positions = tree.inputs.boolean("Pin Positions", False)
        angle = tree.inputs.float("Angle", 0.5983986, subtype="ANGLE")
        distance = tree.inputs.float(
            "Distance", 0.34, min_value=-10_000.0, max_value=10_000.0
        )
        with tree.inputs.panel(
            "Solver",
            description="Solver settings affecting speed and quality.",
            default_closed=True,
        ):
            substeps = tree.inputs.integer(
                "Substeps",
                10,
                description="Number of simulation steps per frame.",
                min_value=1,
                max_value=10000,
                structure_type="SINGLE",
                force_non_field=True,
            )
            constraint_steps = tree.inputs.integer(
                "Constraint Steps",
                15,
                description="Number of steps done to solve constraints in each substep.",
                min_value=1,
                structure_type="SINGLE",
                force_non_field=True,
            )
            time_scale = tree.inputs.float(
                "Time Scale",
                1.0,
                description="Multiplied to the delta time in each time step.",
                min_value=0.0,
                max_value=10.0,
                structure_type="SINGLE",
                force_non_field=True,
            )
            simulation_to_world = tree.inputs.matrix(
                "Simulation to World",
                description="Custom simulation space (by default, simulation space is world space).",
                structure_type="SINGLE",
                force_non_field=True,
            )
        with tree.inputs.panel(
            "Structure",
            description="Settings for structural constraints.",
            default_closed=True,
        ):
            mass = tree.inputs.float(
                "Mass",
                0.01,
                description="Mass assigned to each point in the simulation.",
                min_value=0.001,
                structure_type="FIELD",
                subtype="MASS",
            )
            friction = tree.inputs.float(
                "Friction",
                0.1,
                description="Friction coefficient for each vertex.",
                min_value=0.0,
                structure_type="FIELD",
            )
            _stretchiness = tree.inputs.float(
                "Stretchiness",
                0.0,
                description="Higher values make the hair more stretchy.",
                min_value=0.0,
                max_value=1.0,
                structure_type="FIELD",
                subtype="FACTOR",
            )
            _bendiness = tree.inputs.float(
                "Bendiness",
                0.5,
                description="Higher values allow the hair to bend more.",
                min_value=0.0,
                max_value=1.0,
                structure_type="FIELD",
                subtype="FACTOR",
            )
        with tree.inputs.panel(
            "Damping", description="Damping settings.", default_closed=True
        ):
            linear = tree.inputs.float(
                "Linear",
                3.0,
                description="Simulation stability increases with more linear damping.",
                min_value=0.0,
                structure_type="FIELD",
            )
            angular = tree.inputs.float(
                "Angular",
                3.0,
                description="Simulation stability increases with more angular damping.",
                min_value=0.0,
                structure_type="FIELD",
            )
        with tree.inputs.panel(
            "Surface Collision",
            description="Settings for collisions with the surface mesh.",
            default_closed=True,
        ):
            surface_collision = tree.inputs.boolean(
                "Surface Collision",
                False,
                description="Collide with the surface mesh the hair is attached to.",
                structure_type="SINGLE",
                is_panel_toggle=True,
                force_non_field=True,
            )
            deforming = tree.inputs.boolean(
                "Deforming",
                True,
                description="The surface mesh is deforming.",
                structure_type="SINGLE",
                force_non_field=True,
            )
            edge_contacts = tree.inputs.boolean(
                "Edge Contacts",
                False,
                description="Take edge contacts into account for more accurate collision avoidance.",
                structure_type="SINGLE",
                force_non_field=True,
            )
            surface_friction = tree.inputs.float(
                "Surface Friction",
                0.2,
                description="Friction coefficient of the surface.",
                min_value=0.0,
                max_value=1.0,
                structure_type="SINGLE",
                subtype="FACTOR",
                force_non_field=True,
            )
        with tree.inputs.panel(
            "Effectors",
            description="Additional effectors for the simulation.",
            default_closed=True,
        ):
            effectors_collection = tree.inputs.collection(
                "Effectors Collection",
                description="A collection containing effectors for the simulation. This includes e.g. forces, colliders and custom effectors.",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
            curve_tags = tree.inputs.string(
                "Curve Tags",
                "",
                description="Comma-separated list of tags which can be used to apply only a subset of effectors to this cloth.",
                optional_label=True,
                structure_type="SINGLE",
                force_non_field=True,
            )
            effectors = tree.inputs.bundle(
                "Effectors",
                description="Additional effectors for the simulation.",
                structure_type="SINGLE",
                force_non_field=True,
            )
        dna_curves = tree.outputs.geometry(
            "DNA Curves", description="The deformed hair curves."
        )
        residual_error = tree.outputs.float(
            "Residual Error",
            description="Average relative constraint error. 1 is good, bigger is worse.",
        )

        enable_output = g.EnableOutput.float()
        with g.Frame("Setup Rest Data"):
            get_geometry_component = g.GetGeometryComponent(
                geometry=dna_curve, type="Curve"
            )
            group = SetupStructuralRestData(Geometry=get_geometry_component.o.component)
        with g.Frame("Simulation"):
            with g.Frame("Rod Constraints"):
                with g.Frame("Damping"):
                    combine_bundle = g.CombineBundle()
                    combine_bundle.items.bundle(
                        "Damping", Damping(Linear=linear, Angular=angular)
                    )
                group_1 = RodBendTwistConstraint(
                    Compliance=0.0,
                    **{
                        "Custom Bend Rotation": pin_rotations,
                        "Bend Rotation": g.AxisAngleToRotation(angle=angle),
                    },
                )
                combine_bundle_1 = g.CombineBundle()
                combine_bundle_1.items.bundle(
                    "Stretch_Shear",
                    RodStretchShearConstraint(
                        Compliance=0.0, **{"Custom Length": True}, Length=distance
                    ),
                )
                combine_bundle_1.items.bundle("Bend_Twist", group_1)
            with g.Frame("Surface Collision"):
                group_2 = Collider(
                    Deforming=deforming,
                    **{"Edge Contacts": edge_contacts},
                    Friction=surface_friction,
                )
                combine_bundle_2 = g.CombineBundle()
                combine_bundle_2.items.bundle("Surface Collider", group_2.o.collider)
                switch = surface_collision.switch.bundle(true=combine_bundle_2.o.bundle)
            with g.Frame("Effector Collection"):
                group_3 = CollectionEffector(Collection=effectors_collection)
            with g.Frame("Capture Animated"):
                capture = g.CaptureAttribute.point(geometry=group)
                position = capture.items.vector("Position", g.Position())
                rotation = capture.items.rotation(
                    "Rotation", g.NamedAttribute.quaternion("rotation").o.attribute
                )
            with g.Frame("Curve Root Pinning"):
                endpoint_selection = g.EndpointSelection()
                endpoint_selection_1 = g.EndpointSelection(start_size=2, end_size=2)
                _boolean_math = endpoint_selection.o.selection | pin_rotations
                _boolean_math_1 = endpoint_selection.o.selection | pin_positions
                combine_bundle_3 = g.CombineBundle()
                combine_bundle_3.items.bundle(
                    "Pin Position",
                    PinPositions(
                        Selection=endpoint_selection_1, Position=position.output
                    ),
                )
                combine_bundle_3.items.bundle(
                    "Pin Rotation",
                    PinRotation(
                        Selection=endpoint_selection_1, Rotation=rotation.output
                    ),
                )
            with g.Frame("Setup Curve Geometry for Simulation"):
                group_4 = SetFriction(
                    Geometry=SetGeometryTags(
                        Geometry=capture.o.geometry, Tags=curve_tags
                    ),
                    Static=friction,
                    Dynamic=friction,
                )
                group_5 = SetMass(
                    Geometry=group_4,
                    Mass=mass,
                    **{"Moment of Inertia Mode": "Thin Rod"},
                    _named_links=[("Moment of Inertia", True)],
                )
                combine_bundle_4 = g.CombineBundle()
                combine_bundle_4.items.geometry(
                    "DNA", SetGeoUpdaterDeformOnly(Geometry=group_5)
                )
            join_bundle = g.JoinBundle(
                bundle=(
                    combine_bundle_4.o.bundle,
                    combine_bundle_3.o.bundle,
                    combine_bundle_1.o.bundle,
                    combine_bundle.o.bundle,
                    switch,
                    group_3,
                )
            )
            combine_bundle_5 = g.CombineBundle()
            combine_bundle_5.items.bundle("DNA", join_bundle)
            combine_bundle_5.items.bundle("Effectors", effectors)
            group_6 = XPBDSimulation(
                World=combine_bundle_5.o.bundle,
                Substeps=substeps,
                **{
                    "Constraint Steps": constraint_steps,
                    "Simulation to World": simulation_to_world,
                    "Time Scale": time_scale,
                    "Solver Output Path": "SolverData",
                },
            )
            get_bundle_item = g.GetBundleItem.geometry(group_6, "DNA/DNA")
            _get_bundle_item_1 = g.GetBundleItem.single(
                group_6, "SolverData/residual_error"
            )
        with g.Frame("Transform back to object space"):
            group_7 = ConvertSpaceTransform(
                **{
                    "From Space": "Custom Space",
                    "To Space": "Object Space",
                    "To Object": g.SelfObject(),
                    "Custom to World": simulation_to_world,
                }
            )
            transform_geometry = g.TransformGeometry(
                geometry=get_bundle_item.o.item,
                transform=group_7.o.transform,
                mode="Matrix",
            )
        join_geometry = g.JoinGeometry(
            geometry=(transform_geometry, get_geometry_component.o.geometry)
        )

        join_geometry >> dna_curves
        enable_output >> residual_error


class DNAFromCurve(AssetGeometryGroup):
    """
    DNA From Curve

    Parameters
    ----------
    curves : InputGeometry
        Curves
    base_resolution : InputInteger
        Base Resolution
    socket_3 : InputMenu | Literal["Static", "Simulate"]
        Menu
    wind : InputFloat
        Wind
    socket_5 : InputMenu | Literal["Instance", "Realize"]
        Menu

    Inputs
    ------
    i.curves : GeometrySocket
        Curves
    i.base_resolution : IntegerSocket
        Base Resolution
    i.socket_3 : MenuSocket
        Menu
    i.wind : FloatSocket
        Wind
    i.socket_5 : MenuSocket
        Menu

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "DNA From Curve"
    _asset_name = "DNA From Curve"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        curves: GeometrySocket
        """Curves"""
        base_resolution: IntegerSocket
        """Base Resolution"""
        socket_3: MenuSocket
        """Menu"""
        wind: FloatSocket
        """Wind"""
        socket_5: MenuSocket
        """Menu"""

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
        curves: InputGeometry = None,
        base_resolution: InputInteger = 0,
        socket_3: InputMenu | Literal["Static", "Simulate"] = "Static",
        wind: InputFloat = 1.0,
        socket_5: InputMenu | Literal["Instance", "Realize"] = "Instance",
    ):
        super().__init__(
            **{"Curves": curves, "Base Resolution": base_resolution, "Wind": wind},
            _named_links=[("Menu", socket_3), ("Menu", socket_5)],
        )

    def _build_group(self, tree):
        curves = tree.inputs.geometry("Curves")
        base_resolution = tree.inputs.integer(
            "Base Resolution", 0, min_value=0, max_value=4
        )
        menu = tree.inputs.menu("Menu", expanded=True, optional_label=True)
        wind = tree.inputs.float(
            "Wind", 1.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        menu_1 = tree.inputs.menu("Menu", optional_label=True)
        geometry = tree.outputs.geometry("Geometry")

        with g.Frame("Base instances"):
            group = ColorBackbone(
                backbone=(0.7646205, 0.6673602, 0.7553173, 1.0),
                side_chain=ColorResName(c=(0.20380287, 0.8000075, 0.1007411, 1.0)),
            )
            separate_geometry = SetColor(
                atoms=g.CollectionInfo(
                    collection=bpy.data.collections["DNA Bases"], separate_children=True
                ),
                color=group,
            ) >> g.SeparateGeometry.point(selection=IsHydrogen().o.inverted)
            group_1 = SetColor(
                atoms=separate_geometry.o.selection,
                selection=IsSideChain().o.inverted,
                color=(0.3994269, 0.3737359, 0.3297182, 1.0),
            )
        _group_2 = DNASequenceToID(String="ACGTAAAAAATTAAATTTATATATATATATATATTATATATAT")
        integer_math = base_resolution + 1
        with g.Frame("Compute base to instance"):
            integer_math_1 = ResidueID().o.res_id - 30
            switch = g.NamedAttribute.boolean("is_comp").o.attribute.switch.integer(
                integer_math_1, abs(integer_math_1 - 3)
            )
        math_1 = g.Value(math.tau).o.value / g.Value(10.5)
        accumulate_field = g.Mix(
            factor_float=wind, b_float=math_1 / integer_math, clamp_factor=True
        ).o.result_float.point.trailing(g.CurveOfPoint().o.curve_index)
        store_named_attribute = g.SetHandleType(
            curve=g.SetSplineType.bezier(curves)
        ) >> g.StoreNamedAttribute.point.quaternion(
            name="rotation",
            value=g.AxisAngleToRotation(axis=g.CurveTangent(), angle=accumulate_field),
        )
        with g.Frame("Distance per unwound base"):
            math_2 = g.Value(0.63).o.value * integer_math
        vector_math = g.VectorMath.scale(
            g.NoiseTexture(
                w=g.SceneTime().o.seconds, scale=1.02, noise_dimensions="4D"
            ).o.color,
            4.99,
        )
        group_3 = SimulateDNAGuide(
            **{"DNA Curve": store_named_attribute, "Pin Rotations": True},
            Angle=g.Mix(
                factor_float=wind, b_float=math_1 * integer_math, clamp_factor=True
            ).o.result_float,
            Distance=g.Mix(
                factor_float=wind,
                a_float=math_2,
                b_float=g.Value(0.34),
                clamp_factor=True,
            ).o.result_float,
            Bendiness=0.0,
            Linear=10.0,
            Angular=30.0,
            Effectors=g.JoinBundle(bundle=(CustomForce(Force=vector_math.o.vector),)),
        )
        menu_switch = g.MenuSwitch.geometry(
            menu, {"Static": store_named_attribute, "Simulate": group_3}
        )
        with g.Frame("Handles aren't simulated, so we reset their auto positions"):
            set_handle_type = g.SetHandleType(curve=menu_switch)
        _group_4 = CurveVisualize(
            curve=g.JoinGeometry(geometry=(set_handle_type,)),
            handles=True,
            arrow_size=10.0,
        )
        duplicate_elements = g.DuplicateElements.spline(
            g.SubdivideCurve(curve=set_handle_type, cuts=base_resolution), amount=2
        )
        store_named_attribute_1 = (
            duplicate_elements
            >> g.StoreNamedAttribute.spline.boolean(
                name="is_comp", value=duplicate_elements.o.duplicate_index
            )
        )
        with g.Frame("Compute normal from rotation"):
            axis_angle_to_rotation = g.AxisAngleToRotation(
                angle=g.NamedAttribute.boolean("is_comp").o.attribute.switch.float(
                    true=2.23
                )
            )
            rotate_rotation = g.NamedAttribute.quaternion(
                "rotation"
            ).o.attribute.rotate(axis_angle_to_rotation, rotation_space="LOCAL")
            capture = g.CaptureAttribute.point(geometry=store_named_attribute_1)
            rotation = capture.items.rotation("Rotation", rotate_rotation)
            set_curve_normal = capture.o.geometry >> g.SetCurveNormal(
                normal=g.RotateVector(rotation=rotation.output, vector=(0.0, 1.0, 0.0)),
                mode="Free",
            )
        with g.Frame("Offset from guide curve using normal"):
            capture_1 = g.CaptureAttribute.point(geometry=set_curve_normal)
            capture_1.items.vector("Position", g.Position())
            normal = capture_1.items.vector("Normal", g.Normal().o.normal)
            reverse_curve = (
                capture_1.o.geometry
                >> g.SetPosition(offset=normal.output * AngstromToWorld(angstrom=7.5))
                >> g.SetCurveNormal(normal=normal.output * -1.0, mode="Free")
                >> g.ReverseCurve(
                    selection=g.NamedAttribute.boolean("is_comp").o.attribute
                )
            )
            group_5 = OffsetCurve(curve=reverse_curve)
        capture_2 = g.CaptureAttribute.point(geometry=group_5)
        normal_1 = capture_2.items.vector("Normal", g.Normal().o.normal)
        instance_on_points = capture_2.o.geometry >> g.InstanceOnPoints(
            instance=group_1,
            instance_index=switch,
            rotation=CurveRotation().o.rotation.rotate(
                (0.4262094, 0.23090707, 0.0), rotation_space="LOCAL"
            ),
            pick_instance=True,
        )
        with g.Frame("Rotate individual base nucleotides"):
            capture_3 = g.CaptureAttribute.instance(geometry=instance_on_points)
            transform = capture_3.items.matrix(
                "Transform", TransformLocalAxis(axis=normal_1.output, angle=1.1519172)
            )
            transform_point = g.Position().o.position.transform(
                IsSideChain().o.selection.switch.matrix(true=transform.output)
            )
            set_position = (
                capture_3.o.geometry
                >> g.RealizeInstances(depth=1)
                >> g.SetPosition(position=transform_point)
            )
        (
            g.MenuSwitch.geometry(
                menu_1, {"Instance": instance_on_points, "Realize": set_position}
            )
            >> g.JoinGeometry()
            >> geometry
        )

        menu.default_value = "Static"
        menu_1.default_value = "Instance"


ASSET = DNAFromCurve

ASSET_METADATA = {
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}

DATABLOCK_DEPENDENCIES = {
    "collections": ("DNA Bases",),
    "materials": ("MN Ambient Occlusion",),
}
