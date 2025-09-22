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

    # node Set Material
    set_material = node_group.nodes.new("GeometryNodeSetMaterial")
    set_material.name = "Set Material"
    if "MN Default" in bpy.data.materials:
        set_material.inputs["Material"].default_value = bpy.data.materials["MN Default"]

    # Set locations
    group_input.location = (-340.0, 0.0)
    group_output.location = (560.0, 0.0)
    mesh_to_curve.location = (-160.0, 0.0)
    set_curve_radius.location = (20.0, 0.0)
    thickness_attribute.location = (-160.0, -120.0)
    curve_circle.location = (20.0, -140.0)
    curve_to_mesh.location = (200.0, 0.0)
    set_material.location = (380.0, 0.0)

    # Set dimensions
    group_input.width, group_input.height = 140.0, 100.0
    group_output.width, group_output.height = 140.0, 100.0
    mesh_to_curve.width, mesh_to_curve.height = 140.0, 100.0
    set_curve_radius.width, set_curve_radius.height = 140.0, 100.0
    thickness_attribute.width, thickness_attribute.height = 140.0, 100.0
    curve_circle.width, curve_circle.height = 140.0, 100.0
    curve_to_mesh.width, curve_to_mesh.height = 140.0, 100.0
    set_material.width, set_material.height = 140.0, 100.0

    # initialize links
    # set_material.Geometry -> group_output.Geometry
    node_group.links.new(
        set_material.outputs["Geometry"], group_output.inputs["Geometry"]
    )
    # group_input.Geometry -> mesh_to_curve.Mesh
    node_group.links.new(group_input.outputs["Geometry"], mesh_to_curve.inputs["Mesh"])
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
    # curve_circle.Curve -> curve_to_mesh.Profile Curve
    node_group.links.new(
        curve_circle.outputs["Curve"], curve_to_mesh.inputs["Profile Curve"]
    )
    # curve_to_mesh.Mesh -> set_material.Geometry
    node_group.links.new(curve_to_mesh.outputs["Mesh"], set_material.inputs["Geometry"])
    return node_group
