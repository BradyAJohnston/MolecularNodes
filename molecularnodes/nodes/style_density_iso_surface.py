import bpy


# initialize style_density_iso_surface node group
def style_density_iso_surface_node_group():
    style_density_iso_surface = bpy.data.node_groups.new(
        type="GeometryNodeTree", name="Style Density ISO Surface"
    )

    # style_density_iso_surface interface
    # Socket Geometry
    geometry_socket = style_density_iso_surface.interface.new_socket(
        name="Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry"
    )
    geometry_socket.description = "ISO surface geometry output"

    # Socket Volume
    volume_socket = style_density_iso_surface.interface.new_socket(
        name="Volume", in_out="INPUT", socket_type="NodeSocketGeometry"
    )
    volume_socket.description = "Input geometry"

    # Socket Visible
    visible_socket = style_density_iso_surface.interface.new_socket(
        name="Visible", in_out="INPUT", socket_type="NodeSocketBool"
    )
    visible_socket.default_value = True
    visible_socket.description = "Visibility of style"

    # Socket Shade Smooth
    shade_smooth_socket = style_density_iso_surface.interface.new_socket(
        name="Shade Smooth", in_out="INPUT", socket_type="NodeSocketBool"
    )
    shade_smooth_socket.default_value = True
    shade_smooth_socket.description = "Use smooth shading for surface"

    # Socket ISO Value
    iso_value_socket = style_density_iso_surface.interface.new_socket(
        name="ISO Value", in_out="INPUT", socket_type="NodeSocketFloat"
    )
    iso_value_socket.min_value = 0.0
    iso_value_socket.max_value = 1.0
    iso_value_socket.description = "ISO value"

    # Socket Positive Color
    positive_color_socket = style_density_iso_surface.interface.new_socket(
        name="Positive Color", in_out="INPUT", socket_type="NodeSocketColor"
    )
    positive_color_socket.default_value = (0.0, 0.0, 1.0, 1.0)
    positive_color_socket.description = "Color for positive ISO values"

    # Socket Negative Color
    negative_color_socket = style_density_iso_surface.interface.new_socket(
        name="Negative Color", in_out="INPUT", socket_type="NodeSocketColor"
    )
    negative_color_socket.default_value = (1.0, 0.0, 0.0, 1.0)
    negative_color_socket.description = "Color for negative ISO values"

    # Socket Material
    material_socket = style_density_iso_surface.interface.new_socket(
        name="Material", in_out="INPUT", socket_type="NodeSocketMaterial"
    )
    material_socket.description = "Material to use for this surface"

    # Panel Contours
    contours_panel = style_density_iso_surface.interface.new_panel("Contours")
    # Socket Show Contours
    show_contours_socket = style_density_iso_surface.interface.new_socket(
        name="Show Contours",
        in_out="INPUT",
        socket_type="NodeSocketBool",
        parent=contours_panel,
    )
    show_contours_socket.description = "Whether to show surface contours"

    # Socket Only Contours
    only_contours_socket = style_density_iso_surface.interface.new_socket(
        name="Only Contours",
        in_out="INPUT",
        socket_type="NodeSocketBool",
        parent=contours_panel,
    )
    only_contours_socket.description = "Only show contour edges"

    # Socket Contour Thickness
    contour_thickness_socket = style_density_iso_surface.interface.new_socket(
        name="Contour Thickness",
        in_out="INPUT",
        socket_type="NodeSocketFloat",
        parent=contours_panel,
    )
    contour_thickness_socket.default_value = 0.1
    contour_thickness_socket.min_value = 0.0
    contour_thickness_socket.max_value = 1.0
    contour_thickness_socket.description = "Thickness of the contour edges"

    # Socket Contour Color
    contour_color_socket = style_density_iso_surface.interface.new_socket(
        name="Contour Color",
        in_out="INPUT",
        socket_type="NodeSocketColor",
        parent=contours_panel,
    )
    contour_color_socket.default_value = (0.0, 0.0, 0.0, 1.0)
    contour_color_socket.description = "Color of contour edges"

    # Panel Slicing
    slicing_panel = style_density_iso_surface.interface.new_panel("Slicing")
    # Socket Slice Left
    slice_left_socket = style_density_iso_surface.interface.new_socket(
        name="Slice Left",
        in_out="INPUT",
        socket_type="NodeSocketFloat",
        parent=slicing_panel,
    )
    slice_left_socket.min_value = 0.0
    slice_left_socket.max_value = 100.0
    slice_left_socket.description = "Slice from left (along X axis)"

    # Socket Slice Right
    slice_right_socket = style_density_iso_surface.interface.new_socket(
        name="Slice Right",
        in_out="INPUT",
        socket_type="NodeSocketFloat",
        parent=slicing_panel,
    )
    slice_right_socket.min_value = 0.0
    slice_right_socket.max_value = 100.0
    slice_right_socket.description = "Slice from right (along X axis)"

    # Socket Slice Front
    slice_front_socket = style_density_iso_surface.interface.new_socket(
        name="Slice Front",
        in_out="INPUT",
        socket_type="NodeSocketFloat",
        parent=slicing_panel,
    )
    slice_front_socket.min_value = 0.0
    slice_front_socket.max_value = 100.0
    slice_front_socket.description = "Slice from front (along Y axis)"

    # Socket Slice Back
    slice_back_socket = style_density_iso_surface.interface.new_socket(
        name="Slice Back",
        in_out="INPUT",
        socket_type="NodeSocketFloat",
        parent=slicing_panel,
    )
    slice_back_socket.min_value = 0.0
    slice_back_socket.max_value = 100.0
    slice_back_socket.description = "Slice fom back (along Y axis)"

    # Socket Slice Top
    slice_top_socket = style_density_iso_surface.interface.new_socket(
        name="Slice Top",
        in_out="INPUT",
        socket_type="NodeSocketFloat",
        parent=slicing_panel,
    )
    slice_top_socket.min_value = 0.0
    slice_top_socket.max_value = 100.0
    slice_top_socket.description = "Slice from top (along Z axis)"

    # Socket Slice Bottom
    slice_bottom_socket = style_density_iso_surface.interface.new_socket(
        name="Slice Bottom",
        in_out="INPUT",
        socket_type="NodeSocketFloat",
        parent=slicing_panel,
    )
    slice_bottom_socket.min_value = 0.0
    slice_bottom_socket.max_value = 100.0
    slice_bottom_socket.description = "Slice from bottom (along Z axis)"

    # initialize style_density_iso_surface nodes
    # node Group Output
    group_output = style_density_iso_surface.nodes.new("NodeGroupOutput")
    group_output.name = "Group Output"

    # node Group Input
    group_input = style_density_iso_surface.nodes.new("NodeGroupInput")
    group_input.name = "Group Input"

    # node Mesh to Curve
    mesh_to_curve = style_density_iso_surface.nodes.new("GeometryNodeMeshToCurve")
    mesh_to_curve.name = "Mesh to Curve"

    # node Join Geometry Final
    join_geometry_final = style_density_iso_surface.nodes.new(
        "GeometryNodeJoinGeometry"
    )
    join_geometry_final.name = "Join Geometry Final"

    # node Volume to Mesh Negative
    volume_to_mesh_negative = style_density_iso_surface.nodes.new(
        "GeometryNodeVolumeToMesh"
    )
    volume_to_mesh_negative.name = "Volume to Mesh Negative"

    # node Join Geometry Mesh
    join_geometry_mesh = style_density_iso_surface.nodes.new("GeometryNodeJoinGeometry")
    join_geometry_mesh.name = "Join Geometry Mesh"

    # node Compare Y Positive
    compare_y_positive = style_density_iso_surface.nodes.new("FunctionNodeCompare")
    compare_y_positive.name = "Compare Y Positive"
    compare_y_positive.operation = "LESS_EQUAL"

    # node Delete Geometry
    delete_geometry = style_density_iso_surface.nodes.new("GeometryNodeDeleteGeometry")
    delete_geometry.name = "Delete Geometry"

    # node Compare X Positive
    compare_x_positive = style_density_iso_surface.nodes.new("FunctionNodeCompare")
    compare_x_positive.name = "Compare X Positive"
    compare_x_positive.operation = "LESS_EQUAL"

    # node Position
    position = style_density_iso_surface.nodes.new("GeometryNodeInputPosition")
    position.name = "Position"

    # node Separate XYZ
    separate_xyz = style_density_iso_surface.nodes.new("ShaderNodeSeparateXYZ")
    separate_xyz.name = "Separate XYZ"

    # node Boolean Math X
    boolean_math_x = style_density_iso_surface.nodes.new("FunctionNodeBooleanMath")
    boolean_math_x.name = "Boolean Math X"
    boolean_math_x.operation = "OR"

    # node Set Material Positive
    set_material_positive = style_density_iso_surface.nodes.new(
        "GeometryNodeSetMaterial"
    )
    set_material_positive.name = "Set Material Positive"

    # node Set Material Negative
    set_material_negative = style_density_iso_surface.nodes.new(
        "GeometryNodeSetMaterial"
    )
    set_material_negative.name = "Set Material Negative"

    # node Volume to Mesh Positive
    volume_to_mesh_positive = style_density_iso_surface.nodes.new(
        "GeometryNodeVolumeToMesh"
    )
    volume_to_mesh_positive.name = "Volume to Mesh Positive"

    # node Math
    math = style_density_iso_surface.nodes.new("ShaderNodeMath")
    math.name = "Math"
    math.operation = "MULTIPLY"
    # Value_001
    math.inputs[1].default_value = -1.0

    # node Set Shade Smooth
    set_shade_smooth = style_density_iso_surface.nodes.new("GeometryNodeSetShadeSmooth")
    set_shade_smooth.name = "Set Shade Smooth"

    # node Compare Z Positive
    compare_z_positive = style_density_iso_surface.nodes.new("FunctionNodeCompare")
    compare_z_positive.name = "Compare Z Positive"
    compare_z_positive.operation = "LESS_EQUAL"

    # node Boolean Math XY
    boolean_math_xy = style_density_iso_surface.nodes.new("FunctionNodeBooleanMath")
    boolean_math_xy.name = "Boolean Math XY"
    boolean_math_xy.operation = "OR"

    # node Map Range X Positive
    map_range_x_positive = style_density_iso_surface.nodes.new("ShaderNodeMapRange")
    map_range_x_positive.name = "Map Range X Positive"
    # From Max
    map_range_x_positive.inputs[2].default_value = 100.0

    # node Map Range Y Positive
    map_range_y_positive = style_density_iso_surface.nodes.new("ShaderNodeMapRange")
    map_range_y_positive.name = "Map Range Y Positive"
    # From Max
    map_range_y_positive.inputs[2].default_value = 100.0

    # node Map Range Z Positive
    map_range_z_positive = style_density_iso_surface.nodes.new("ShaderNodeMapRange")
    map_range_z_positive.name = "Map Range Z Positive"
    # From Max
    map_range_z_positive.inputs[2].default_value = 100.0

    # node Map Range X Negative
    map_range_x_negative = style_density_iso_surface.nodes.new("ShaderNodeMapRange")
    map_range_x_negative.name = "Map Range X Negative"
    # From Max
    map_range_x_negative.inputs[2].default_value = 100.0

    # node Compare X Negative
    compare_x_negative = style_density_iso_surface.nodes.new("FunctionNodeCompare")
    compare_x_negative.name = "Compare X Negative"
    compare_x_negative.operation = "GREATER_EQUAL"

    # node Boolean Math XYZ
    boolean_math_xyz = style_density_iso_surface.nodes.new("FunctionNodeBooleanMath")
    boolean_math_xyz.name = "Boolean Math XYZ"
    boolean_math_xyz.operation = "OR"

    # node X Min
    x_min = style_density_iso_surface.nodes.new("ShaderNodeValue")
    x_min.name = "X Min"

    # node X Max
    x_max = style_density_iso_surface.nodes.new("ShaderNodeValue")
    x_max.name = "X Max"

    x_max.outputs["Value"].default_value = 1.0
    # node Frame X
    frame_x = style_density_iso_surface.nodes.new("NodeFrame")
    frame_x.name = "Frame X"

    # node Map Range Y Negative
    map_range_y_negative = style_density_iso_surface.nodes.new("ShaderNodeMapRange")
    map_range_y_negative.name = "Map Range Y Negative"
    # From Max
    map_range_y_negative.inputs[2].default_value = 100.0

    # node Y Min
    y_min = style_density_iso_surface.nodes.new("ShaderNodeValue")
    y_min.name = "Y Min"

    # node Y Max
    y_max = style_density_iso_surface.nodes.new("ShaderNodeValue")
    y_max.name = "Y Max"

    y_max.outputs["Value"].default_value = 1.0
    # node Compare Y Negative
    compare_y_negative = style_density_iso_surface.nodes.new("FunctionNodeCompare")
    compare_y_negative.name = "Compare Y Negative"
    compare_y_negative.operation = "GREATER_EQUAL"

    # node Boolean Math Y
    boolean_math_y = style_density_iso_surface.nodes.new("FunctionNodeBooleanMath")
    boolean_math_y.name = "Boolean Math Y"
    boolean_math_y.operation = "OR"

    # node Frame Y
    frame_y = style_density_iso_surface.nodes.new("NodeFrame")
    frame_y.name = "Frame Y"

    # node Z Max
    z_max = style_density_iso_surface.nodes.new("ShaderNodeValue")
    z_max.name = "Z Max"

    # node Z Min
    z_min = style_density_iso_surface.nodes.new("ShaderNodeValue")
    z_min.name = "Z Min"

    z_min.outputs["Value"].default_value = 1.0
    # node Map Range Z Negative
    map_range_z_negative = style_density_iso_surface.nodes.new("ShaderNodeMapRange")
    map_range_z_negative.name = "Map Range Z Negative"
    # From Max
    map_range_z_negative.inputs[2].default_value = 100.0

    # node Compare Z Negative
    compare_z_negative = style_density_iso_surface.nodes.new("FunctionNodeCompare")
    compare_z_negative.name = "Compare Z Negative"
    compare_z_negative.operation = "GREATER_EQUAL"

    # node Boolean Math Z
    boolean_math_z = style_density_iso_surface.nodes.new("FunctionNodeBooleanMath")
    boolean_math_z.name = "Boolean Math Z"
    boolean_math_z.operation = "OR"

    # node Frame Z
    frame_z = style_density_iso_surface.nodes.new("NodeFrame")
    frame_z.name = "Frame Z"

    # node Store Named Attribute Positive
    store_named_attribute_positive = style_density_iso_surface.nodes.new(
        "GeometryNodeStoreNamedAttribute"
    )
    store_named_attribute_positive.name = "Store Named Attribute Positive"
    store_named_attribute_positive.data_type = "FLOAT_COLOR"
    # Name
    store_named_attribute_positive.inputs["Name"].default_value = "Color"

    # node Store Named Attribute Negative
    store_named_attribute_negative = style_density_iso_surface.nodes.new(
        "GeometryNodeStoreNamedAttribute"
    )
    store_named_attribute_negative.name = "Store Named Attribute Negative"
    store_named_attribute_negative.data_type = "FLOAT_COLOR"
    # Name
    store_named_attribute_negative.inputs["Name"].default_value = "Color"

    # node Curve to Mesh
    curve_to_mesh = style_density_iso_surface.nodes.new("GeometryNodeCurveToMesh")
    curve_to_mesh.name = "Curve to Mesh"

    # node Quadrilateral
    quadrilateral = style_density_iso_surface.nodes.new(
        "GeometryNodeCurvePrimitiveQuadrilateral"
    )
    quadrilateral.name = "Quadrilateral"

    # node Store Named Attribute Contours
    store_named_attribute_contours = style_density_iso_surface.nodes.new(
        "GeometryNodeStoreNamedAttribute"
    )
    store_named_attribute_contours.name = "Store Named Attribute Contours"
    store_named_attribute_contours.data_type = "FLOAT_COLOR"
    store_named_attribute_contours.domain = "EDGE"
    # Name
    store_named_attribute_contours.inputs["Name"].default_value = "Color"

    # node Set Material Contours
    set_material_contours = style_density_iso_surface.nodes.new(
        "GeometryNodeSetMaterial"
    )
    set_material_contours.name = "Set Material Contours"

    # node Scale Down Contour Thickness
    scale_down_contour_thickness = style_density_iso_surface.nodes.new("ShaderNodeMath")
    scale_down_contour_thickness.name = "Scale Down Contour Thickness"
    scale_down_contour_thickness.operation = "MULTIPLY"
    # Value_001
    scale_down_contour_thickness.inputs[1].default_value = 0.001

    # node Only Contours Switch
    only_contours_switch = style_density_iso_surface.nodes.new("GeometryNodeSwitch")
    only_contours_switch.name = "Only Contours Switch"

    # node Visibility
    visibility = style_density_iso_surface.nodes.new("GeometryNodeSwitch")
    visibility.name = "Visibility"

    # Set parents
    compare_y_positive.parent = frame_y
    compare_x_positive.parent = frame_x
    boolean_math_x.parent = frame_x
    compare_z_positive.parent = frame_z
    map_range_x_positive.parent = frame_x
    map_range_y_positive.parent = frame_y
    map_range_z_positive.parent = frame_z
    map_range_x_negative.parent = frame_x
    compare_x_negative.parent = frame_x
    x_min.parent = frame_x
    x_max.parent = frame_x
    map_range_y_negative.parent = frame_y
    y_min.parent = frame_y
    y_max.parent = frame_y
    compare_y_negative.parent = frame_y
    boolean_math_y.parent = frame_y
    z_max.parent = frame_z
    z_min.parent = frame_z
    map_range_z_negative.parent = frame_z
    compare_z_negative.parent = frame_z
    boolean_math_z.parent = frame_z

    # Set locations
    group_output.location = (2200.0, 40.0)
    group_input.location = (-1100.0, 0.0)
    mesh_to_curve.location = (1180.0, -80.0)
    join_geometry_final.location = (2000.0, 40.0)
    volume_to_mesh_negative.location = (-400.0, 20.0)
    join_geometry_mesh.location = (380.0, 160.0)
    compare_y_positive.location = (550.0, -30.0)
    delete_geometry.location = (1000.0, 20.0)
    compare_x_positive.location = (550.0, -30.0)
    position.location = (-940.0, -700.0)
    separate_xyz.location = (-760.0, -680.0)
    boolean_math_x.location = (730.0, -90.0)
    set_material_positive.location = (-60.0, 220.0)
    set_material_negative.location = (0.0, 40.0)
    volume_to_mesh_positive.location = (-460.0, 200.0)
    math.location = (-580.0, 20.0)
    set_shade_smooth.location = (580.0, 200.0)
    compare_z_positive.location = (530.0, -30.0)
    boolean_math_xy.location = (520.0, -360.0)
    map_range_x_positive.location = (210.0, -70.0)
    map_range_y_positive.location = (210.0, -70.0)
    map_range_z_positive.location = (190.0, -50.0)
    map_range_x_negative.location = (370.0, -70.0)
    compare_x_negative.location = (550.0, -190.0)
    boolean_math_xyz.location = (720.0, -440.0)
    x_min.location = (30.0, -110.0)
    x_max.location = (30.0, -210.0)
    frame_x.location = (-450.0, -190.0)
    map_range_y_negative.location = (370.0, -70.0)
    y_min.location = (30.0, -90.0)
    y_max.location = (30.0, -190.0)
    compare_y_negative.location = (550.0, -190.0)
    boolean_math_y.location = (730.0, -70.0)
    frame_y.location = (-450.0, -590.0)
    z_max.location = (30.0, -170.0)
    z_min.location = (30.0, -70.0)
    map_range_z_negative.location = (350.0, -50.0)
    compare_z_negative.location = (530.0, -190.0)
    boolean_math_z.location = (710.0, -70.0)
    frame_z.location = (-430.0, -990.0)
    store_named_attribute_positive.location = (-240.0, 260.0)
    store_named_attribute_negative.location = (-180.0, 40.0)
    curve_to_mesh.location = (1360.0, -80.0)
    quadrilateral.location = (1180.0, -200.0)
    store_named_attribute_contours.location = (1540.0, -80.0)
    set_material_contours.location = (1720.0, -40.0)
    scale_down_contour_thickness.location = (1020.0, -200.0)
    only_contours_switch.location = (1360.0, 80.0)
    visibility.location = (-840.0, 120.0)

    # Set dimensions
    frame_x.width, frame_x.height = 900.0, 368.0
    frame_y.width, frame_y.height = 900.0, 368.0
    frame_z.width, frame_z.height = 880.0, 368.0

    # initialize style_density_iso_surface links
    # delete_geometry.Geometry -> mesh_to_curve.Mesh
    style_density_iso_surface.links.new(
        delete_geometry.outputs["Geometry"], mesh_to_curve.inputs["Mesh"]
    )
    # separate_xyz.X -> compare_x_positive.A
    style_density_iso_surface.links.new(
        separate_xyz.outputs["X"], compare_x_positive.inputs[0]
    )
    # math.Value -> volume_to_mesh_negative.Threshold
    style_density_iso_surface.links.new(
        math.outputs["Value"], volume_to_mesh_negative.inputs["Threshold"]
    )
    # store_named_attribute_positive.Geometry -> set_material_positive.Geometry
    style_density_iso_surface.links.new(
        store_named_attribute_positive.outputs["Geometry"],
        set_material_positive.inputs["Geometry"],
    )
    # join_geometry_mesh.Geometry -> set_shade_smooth.Geometry
    style_density_iso_surface.links.new(
        join_geometry_mesh.outputs["Geometry"], set_shade_smooth.inputs["Geometry"]
    )
    # separate_xyz.Y -> compare_y_positive.A
    style_density_iso_surface.links.new(
        separate_xyz.outputs["Y"], compare_y_positive.inputs[0]
    )
    # store_named_attribute_negative.Geometry -> set_material_negative.Geometry
    style_density_iso_surface.links.new(
        store_named_attribute_negative.outputs["Geometry"],
        set_material_negative.inputs["Geometry"],
    )
    # set_shade_smooth.Geometry -> delete_geometry.Geometry
    style_density_iso_surface.links.new(
        set_shade_smooth.outputs["Geometry"], delete_geometry.inputs["Geometry"]
    )
    # position.Position -> separate_xyz.Vector
    style_density_iso_surface.links.new(
        position.outputs["Position"], separate_xyz.inputs["Vector"]
    )
    # set_material_negative.Geometry -> join_geometry_mesh.Geometry
    style_density_iso_surface.links.new(
        set_material_negative.outputs["Geometry"], join_geometry_mesh.inputs["Geometry"]
    )
    # compare_x_positive.Result -> boolean_math_x.Boolean
    style_density_iso_surface.links.new(
        compare_x_positive.outputs["Result"], boolean_math_x.inputs[0]
    )
    # group_input.ISO Value -> volume_to_mesh_positive.Threshold
    style_density_iso_surface.links.new(
        group_input.outputs["ISO Value"], volume_to_mesh_positive.inputs["Threshold"]
    )
    # group_input.ISO Value -> math.Value
    style_density_iso_surface.links.new(
        group_input.outputs["ISO Value"], math.inputs[0]
    )
    # group_input.Show Contours -> mesh_to_curve.Selection
    style_density_iso_surface.links.new(
        group_input.outputs["Show Contours"], mesh_to_curve.inputs["Selection"]
    )
    # visibility.Output -> volume_to_mesh_positive.Volume
    style_density_iso_surface.links.new(
        visibility.outputs["Output"], volume_to_mesh_positive.inputs["Volume"]
    )
    # join_geometry_final.Geometry -> group_output.Geometry
    style_density_iso_surface.links.new(
        join_geometry_final.outputs["Geometry"], group_output.inputs["Geometry"]
    )
    # boolean_math_x.Boolean -> boolean_math_xy.Boolean
    style_density_iso_surface.links.new(
        boolean_math_x.outputs["Boolean"], boolean_math_xy.inputs[0]
    )
    # group_input.Slice Left -> map_range_x_positive.Value
    style_density_iso_surface.links.new(
        group_input.outputs["Slice Left"], map_range_x_positive.inputs["Value"]
    )
    # group_input.Slice Front -> map_range_y_positive.Value
    style_density_iso_surface.links.new(
        group_input.outputs["Slice Front"], map_range_y_positive.inputs["Value"]
    )
    # group_input.Slice Top -> map_range_z_negative.Value
    style_density_iso_surface.links.new(
        group_input.outputs["Slice Top"], map_range_z_negative.inputs["Value"]
    )
    # map_range_y_positive.Result -> compare_y_positive.B
    style_density_iso_surface.links.new(
        map_range_y_positive.outputs["Result"], compare_y_positive.inputs[1]
    )
    # separate_xyz.Z -> compare_z_positive.A
    style_density_iso_surface.links.new(
        separate_xyz.outputs["Z"], compare_z_positive.inputs[0]
    )
    # map_range_z_positive.Result -> compare_z_positive.B
    style_density_iso_surface.links.new(
        map_range_z_positive.outputs["Result"], compare_z_positive.inputs[1]
    )
    # map_range_x_positive.Result -> compare_x_positive.B
    style_density_iso_surface.links.new(
        map_range_x_positive.outputs["Result"], compare_x_positive.inputs[1]
    )
    # group_input.Slice Right -> map_range_x_negative.Value
    style_density_iso_surface.links.new(
        group_input.outputs["Slice Right"], map_range_x_negative.inputs["Value"]
    )
    # separate_xyz.X -> compare_x_negative.A
    style_density_iso_surface.links.new(
        separate_xyz.outputs["X"], compare_x_negative.inputs[0]
    )
    # map_range_x_negative.Result -> compare_x_negative.B
    style_density_iso_surface.links.new(
        map_range_x_negative.outputs["Result"], compare_x_negative.inputs[1]
    )
    # compare_x_negative.Result -> boolean_math_x.Boolean
    style_density_iso_surface.links.new(
        compare_x_negative.outputs["Result"], boolean_math_x.inputs[1]
    )
    # boolean_math_z.Boolean -> boolean_math_xyz.Boolean
    style_density_iso_surface.links.new(
        boolean_math_z.outputs["Boolean"], boolean_math_xyz.inputs[1]
    )
    # x_min.Value -> map_range_x_positive.To Min
    style_density_iso_surface.links.new(
        x_min.outputs["Value"], map_range_x_positive.inputs[3]
    )
    # x_max.Value -> map_range_x_positive.To Max
    style_density_iso_surface.links.new(
        x_max.outputs["Value"], map_range_x_positive.inputs[4]
    )
    # x_min.Value -> map_range_x_negative.To Max
    style_density_iso_surface.links.new(
        x_min.outputs["Value"], map_range_x_negative.inputs[4]
    )
    # x_max.Value -> map_range_x_negative.To Min
    style_density_iso_surface.links.new(
        x_max.outputs["Value"], map_range_x_negative.inputs[3]
    )
    # y_min.Value -> map_range_y_positive.To Min
    style_density_iso_surface.links.new(
        y_min.outputs["Value"], map_range_y_positive.inputs[3]
    )
    # y_max.Value -> map_range_y_positive.To Max
    style_density_iso_surface.links.new(
        y_max.outputs["Value"], map_range_y_positive.inputs[4]
    )
    # y_min.Value -> map_range_y_negative.To Max
    style_density_iso_surface.links.new(
        y_min.outputs["Value"], map_range_y_negative.inputs[4]
    )
    # y_max.Value -> map_range_y_negative.To Min
    style_density_iso_surface.links.new(
        y_max.outputs["Value"], map_range_y_negative.inputs[3]
    )
    # separate_xyz.Y -> compare_y_negative.A
    style_density_iso_surface.links.new(
        separate_xyz.outputs["Y"], compare_y_negative.inputs[0]
    )
    # map_range_y_negative.Result -> compare_y_negative.B
    style_density_iso_surface.links.new(
        map_range_y_negative.outputs["Result"], compare_y_negative.inputs[1]
    )
    # compare_y_positive.Result -> boolean_math_y.Boolean
    style_density_iso_surface.links.new(
        compare_y_positive.outputs["Result"], boolean_math_y.inputs[0]
    )
    # compare_y_negative.Result -> boolean_math_y.Boolean
    style_density_iso_surface.links.new(
        compare_y_negative.outputs["Result"], boolean_math_y.inputs[1]
    )
    # group_input.Slice Back -> map_range_y_negative.Value
    style_density_iso_surface.links.new(
        group_input.outputs["Slice Back"], map_range_y_negative.inputs["Value"]
    )
    # z_min.Value -> map_range_z_positive.To Min
    style_density_iso_surface.links.new(
        z_min.outputs["Value"], map_range_z_positive.inputs[3]
    )
    # z_max.Value -> map_range_z_positive.To Max
    style_density_iso_surface.links.new(
        z_max.outputs["Value"], map_range_z_positive.inputs[4]
    )
    # z_min.Value -> map_range_z_negative.To Max
    style_density_iso_surface.links.new(
        z_min.outputs["Value"], map_range_z_negative.inputs[4]
    )
    # z_max.Value -> map_range_z_negative.To Min
    style_density_iso_surface.links.new(
        z_max.outputs["Value"], map_range_z_negative.inputs[3]
    )
    # group_input.Slice Bottom -> map_range_z_positive.Value
    style_density_iso_surface.links.new(
        group_input.outputs["Slice Bottom"], map_range_z_positive.inputs["Value"]
    )
    # map_range_z_negative.Result -> compare_z_negative.B
    style_density_iso_surface.links.new(
        map_range_z_negative.outputs["Result"], compare_z_negative.inputs[1]
    )
    # separate_xyz.Z -> compare_z_negative.A
    style_density_iso_surface.links.new(
        separate_xyz.outputs["Z"], compare_z_negative.inputs[0]
    )
    # compare_z_positive.Result -> boolean_math_z.Boolean
    style_density_iso_surface.links.new(
        compare_z_positive.outputs["Result"], boolean_math_z.inputs[0]
    )
    # compare_z_negative.Result -> boolean_math_z.Boolean
    style_density_iso_surface.links.new(
        compare_z_negative.outputs["Result"], boolean_math_z.inputs[1]
    )
    # boolean_math_y.Boolean -> boolean_math_xy.Boolean
    style_density_iso_surface.links.new(
        boolean_math_y.outputs["Boolean"], boolean_math_xy.inputs[1]
    )
    # boolean_math_xy.Boolean -> boolean_math_xyz.Boolean
    style_density_iso_surface.links.new(
        boolean_math_xy.outputs["Boolean"], boolean_math_xyz.inputs[0]
    )
    # boolean_math_xyz.Boolean -> delete_geometry.Selection
    style_density_iso_surface.links.new(
        boolean_math_xyz.outputs["Boolean"], delete_geometry.inputs["Selection"]
    )
    # group_input.Shade Smooth -> set_shade_smooth.Shade Smooth
    style_density_iso_surface.links.new(
        group_input.outputs["Shade Smooth"], set_shade_smooth.inputs["Shade Smooth"]
    )
    # volume_to_mesh_positive.Mesh -> store_named_attribute_positive.Geometry
    style_density_iso_surface.links.new(
        volume_to_mesh_positive.outputs["Mesh"],
        store_named_attribute_positive.inputs["Geometry"],
    )
    # volume_to_mesh_negative.Mesh -> store_named_attribute_negative.Geometry
    style_density_iso_surface.links.new(
        volume_to_mesh_negative.outputs["Mesh"],
        store_named_attribute_negative.inputs["Geometry"],
    )
    # group_input.Positive Color -> store_named_attribute_positive.Value
    style_density_iso_surface.links.new(
        group_input.outputs["Positive Color"],
        store_named_attribute_positive.inputs["Value"],
    )
    # group_input.Negative Color -> store_named_attribute_negative.Value
    style_density_iso_surface.links.new(
        group_input.outputs["Negative Color"],
        store_named_attribute_negative.inputs["Value"],
    )
    # group_input.Material -> set_material_positive.Material
    style_density_iso_surface.links.new(
        group_input.outputs["Material"], set_material_positive.inputs["Material"]
    )
    # group_input.Material -> set_material_negative.Material
    style_density_iso_surface.links.new(
        group_input.outputs["Material"], set_material_negative.inputs["Material"]
    )
    # set_material_contours.Geometry -> join_geometry_final.Geometry
    style_density_iso_surface.links.new(
        set_material_contours.outputs["Geometry"],
        join_geometry_final.inputs["Geometry"],
    )
    # mesh_to_curve.Curve -> curve_to_mesh.Curve
    style_density_iso_surface.links.new(
        mesh_to_curve.outputs["Curve"], curve_to_mesh.inputs["Curve"]
    )
    # quadrilateral.Curve -> curve_to_mesh.Profile Curve
    style_density_iso_surface.links.new(
        quadrilateral.outputs["Curve"], curve_to_mesh.inputs["Profile Curve"]
    )
    # curve_to_mesh.Mesh -> store_named_attribute_contours.Geometry
    style_density_iso_surface.links.new(
        curve_to_mesh.outputs["Mesh"], store_named_attribute_contours.inputs["Geometry"]
    )
    # store_named_attribute_contours.Geometry -> set_material_contours.Geometry
    style_density_iso_surface.links.new(
        store_named_attribute_contours.outputs["Geometry"],
        set_material_contours.inputs["Geometry"],
    )
    # group_input.Material -> set_material_contours.Material
    style_density_iso_surface.links.new(
        group_input.outputs["Material"], set_material_contours.inputs["Material"]
    )
    # scale_down_contour_thickness.Value -> quadrilateral.Width
    style_density_iso_surface.links.new(
        scale_down_contour_thickness.outputs["Value"], quadrilateral.inputs["Width"]
    )
    # scale_down_contour_thickness.Value -> quadrilateral.Height
    style_density_iso_surface.links.new(
        scale_down_contour_thickness.outputs["Value"], quadrilateral.inputs["Height"]
    )
    # group_input.Contour Thickness -> scale_down_contour_thickness.Value
    style_density_iso_surface.links.new(
        group_input.outputs["Contour Thickness"], scale_down_contour_thickness.inputs[0]
    )
    # delete_geometry.Geometry -> only_contours_switch.False
    style_density_iso_surface.links.new(
        delete_geometry.outputs["Geometry"], only_contours_switch.inputs["False"]
    )
    # group_input.Only Contours -> only_contours_switch.Switch
    style_density_iso_surface.links.new(
        group_input.outputs["Only Contours"], only_contours_switch.inputs["Switch"]
    )
    # group_input.Contour Color -> store_named_attribute_contours.Value
    style_density_iso_surface.links.new(
        group_input.outputs["Contour Color"],
        store_named_attribute_contours.inputs["Value"],
    )
    # group_input.Volume -> visibility.True
    style_density_iso_surface.links.new(
        group_input.outputs["Volume"], visibility.inputs["True"]
    )
    # visibility.Output -> volume_to_mesh_negative.Volume
    style_density_iso_surface.links.new(
        visibility.outputs["Output"], volume_to_mesh_negative.inputs["Volume"]
    )
    # group_input.Visible -> visibility.Switch
    style_density_iso_surface.links.new(
        group_input.outputs["Visible"], visibility.inputs["Switch"]
    )
    # set_material_positive.Geometry -> join_geometry_mesh.Geometry
    style_density_iso_surface.links.new(
        set_material_positive.outputs["Geometry"], join_geometry_mesh.inputs["Geometry"]
    )
    # only_contours_switch.Output -> join_geometry_final.Geometry
    style_density_iso_surface.links.new(
        only_contours_switch.outputs["Output"], join_geometry_final.inputs["Geometry"]
    )
    return style_density_iso_surface
