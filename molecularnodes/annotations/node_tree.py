import bpy


def annotations_node_tree():
    node_group = bpy.data.node_groups.new(type="GeometryNodeTree", name="Annotations")
    node_group.is_modifier = True

    # Output Socket Geometry
    node_group.interface.new_socket(
        name="Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry"
    )

    # Input Socket Geometry
    node_group.interface.new_socket(
        name="Geometry", in_out="INPUT", socket_type="NodeSocketGeometry"
    )

    # initialize nodes
    # node Group Input
    group_input = node_group.nodes.new("NodeGroupInput")
    group_input.name = "Group Input"

    # node Group Output
    group_output = node_group.nodes.new("NodeGroupOutput")
    group_output.name = "Group Output"

    # node Mesh to Curve
    mesh_to_curve = node_group.nodes.new("GeometryNodeMeshToCurve")
    mesh_to_curve.name = "Mesh to Curve"

    # node Set Curve Radius
    set_curve_radius = node_group.nodes.new("GeometryNodeSetCurveRadius")
    set_curve_radius.name = "Set Curve Radius"

    # node Is Line Attribute
    is_line_attribute = node_group.nodes.new("GeometryNodeInputNamedAttribute")
    is_line_attribute.name = "Is Line Attribute"
    is_line_attribute.data_type = "BOOLEAN"
    is_line_attribute.inputs["Name"].default_value = "is_line"

    # node Thickness Attribute
    thickness_attribute = node_group.nodes.new("GeometryNodeInputNamedAttribute")
    thickness_attribute.name = "Thickness Attribute"
    thickness_attribute.inputs["Name"].default_value = "thickness"

    # node Curve Circle
    curve_circle = node_group.nodes.new("GeometryNodeCurvePrimitiveCircle")
    curve_circle.name = "Curve Circle"
    curve_circle.inputs["Resolution"].default_value = 32
    curve_circle.inputs["Radius"].default_value = 0.001

    # node Curve to Mesh
    curve_to_mesh = node_group.nodes.new("GeometryNodeCurveToMesh")
    curve_to_mesh.name = "Curve to Mesh"
    curve_to_mesh.inputs["Fill Caps"].default_value = True

    # node Join Geometry
    join_geometry = node_group.nodes.new("GeometryNodeJoinGeometry")
    join_geometry.name = "Join Geometry"

    # node Set Material Index
    set_material_index = node_group.nodes.new("GeometryNodeSetMaterialIndex")
    set_material_index.name = "Set Material Index"

    # node Material Slot Attribute
    material_slot_attribute = node_group.nodes.new("GeometryNodeInputNamedAttribute")
    material_slot_attribute.name = "Material Slot Attribute"
    material_slot_attribute.data_type = "INT"
    # Name
    material_slot_attribute.inputs["Name"].default_value = "material_slot_index"

    # node Shade Smooth Attribute
    shade_smooth_attribute = node_group.nodes.new("GeometryNodeInputNamedAttribute")
    shade_smooth_attribute.name = "Shade Smooth Attribute"
    shade_smooth_attribute.data_type = "BOOLEAN"
    # Name
    shade_smooth_attribute.inputs["Name"].default_value = "shade_smooth"

    # node Set Shade Smooth
    set_shade_smooth = node_group.nodes.new("GeometryNodeSetShadeSmooth")
    set_shade_smooth.name = "Set Shade Smooth"

    # Set locations
    group_input.location = (-520.0, 0.0)
    group_output.location = (920.0, 0.0)
    mesh_to_curve.location = (-160.0, 0.0)
    set_curve_radius.location = (20.0, 0.0)
    is_line_attribute.location = (-340.0, -120.0)
    thickness_attribute.location = (-160.0, -120.0)
    curve_circle.location = (20.0, -140.0)
    curve_to_mesh.location = (200.0, 0.0)
    join_geometry.location = (380.0, 0.0)
    set_material_index.location = (560.0, 0.0)
    material_slot_attribute.location = (380.0, -120.0)
    shade_smooth_attribute.location = (560.0, -140.0)
    set_shade_smooth.location = (740.0, 0.0)

    # initialize links
    # is_line_attribute.Attribute -> set_shade_smooth.Selection
    node_group.links.new(
        shade_smooth_attribute.outputs["Attribute"],
        set_shade_smooth.inputs["Shade Smooth"],
    )
    # set_material_index.Geometry -> set_shade_smooth.Geometry
    node_group.links.new(
        set_material_index.outputs["Geometry"], set_shade_smooth.inputs["Geometry"]
    )
    # set_shade_smooth.Geometry -> group_output.Geometry
    node_group.links.new(
        set_shade_smooth.outputs["Geometry"], group_output.inputs["Geometry"]
    )
    # join_geometry.Geometry -> set_material_index.Geometry
    node_group.links.new(
        join_geometry.outputs["Geometry"], set_material_index.inputs["Geometry"]
    )
    # material_slot_attribute.Attribute -> set_material_index.Material Index
    node_group.links.new(
        material_slot_attribute.outputs["Attribute"],
        set_material_index.inputs["Material Index"],
    )
    # group_input.Geometry -> mesh_to_curve.Mesh
    node_group.links.new(group_input.outputs["Geometry"], mesh_to_curve.inputs["Mesh"])
    if bpy.app.version >= (4, 5, 0):
        # From Blender 4.5.x onwards, need to use Scale input of Curve to Mesh
        # remove set_curve_radius node
        node_group.nodes.remove(set_curve_radius)
        # mesh_to_curve.Curve -> curve_to_mesh.Curve
        node_group.links.new(
            mesh_to_curve.outputs["Curve"], curve_to_mesh.inputs["Curve"]
        )
        # thickness_attribute.Attribute -> curve_to_mesh.Scale
        node_group.links.new(
            thickness_attribute.outputs["Attribute"], curve_to_mesh.inputs["Scale"]
        )
    else:
        # mesh_to_curve.Curve -> set_curve_radius.Curve
        node_group.links.new(
            mesh_to_curve.outputs["Curve"], set_curve_radius.inputs["Curve"]
        )
        # thickness_attribute.Attribute -> set_curve_radius.Radius
        node_group.links.new(
            thickness_attribute.outputs["Attribute"], set_curve_radius.inputs["Radius"]
        )
        # set_curve_radius.Curve -> curve_to_mesh.Curve
        node_group.links.new(
            set_curve_radius.outputs["Curve"], curve_to_mesh.inputs["Curve"]
        )
    # is_line_attribute.Attribute -> mesh_to_curve.Selection
    node_group.links.new(
        is_line_attribute.outputs["Attribute"], mesh_to_curve.inputs["Selection"]
    )
    # curve_circle.Curve -> curve_to_mesh.Profile Curve
    node_group.links.new(
        curve_circle.outputs["Curve"], curve_to_mesh.inputs["Profile Curve"]
    )
    # curve_to_mesh.Mesh -> join_geometry.Geometry
    node_group.links.new(
        curve_to_mesh.outputs["Mesh"], join_geometry.inputs["Geometry"]
    )
    # group_input.Geometry -> join_geometry.Geometry
    node_group.links.new(
        group_input.outputs["Geometry"], join_geometry.inputs["Geometry"]
    )
    return node_group
