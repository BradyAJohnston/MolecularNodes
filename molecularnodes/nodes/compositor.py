import bpy


# Default 5.x compositor node tree
def default_5x_compositor_node_tree():
    node_tree = bpy.data.node_groups.new(
        type="CompositorNodeTree", name="Compositor Nodes"
    )
    render_layers = node_tree.nodes.new("CompositorNodeRLayers")
    group_output = node_tree.nodes.new("NodeGroupOutput")
    _ = node_tree.interface.new_socket(
        name="Image", in_out="OUTPUT", socket_type="NodeSocketColor"
    )
    render_layers.location = (100.0, 250.0)
    group_output.location = (600.0, 250.0)
    node_tree.links.new(render_layers.outputs["Image"], group_output.inputs["Image"])
    return node_tree


# MN Compositor node tree
def mn_compositor_node_tree():
    mn_compositor = bpy.data.node_groups.new(
        type="CompositorNodeTree", name="MN Compositor"
    )

    # mn_compositor interface
    # Socket Image Output
    mn_compositor.interface.new_socket(
        name="Image", in_out="OUTPUT", socket_type="NodeSocketColor"
    )
    # Socket Image Input
    mn_compositor.interface.new_socket(
        name="Image", in_out="INPUT", socket_type="NodeSocketColor"
    )

    # initialize mn_compositor nodes
    # node Group Output
    group_output = mn_compositor.nodes.new("NodeGroupOutput")
    group_output.name = "Group Output"

    # node Group Input
    group_input = mn_compositor.nodes.new("NodeGroupInput")
    group_input.name = "Group Input"

    # node Image
    image = mn_compositor.nodes.new("CompositorNodeImage")
    image.name = "Image"
    if "mn_annotations" in bpy.data.images:
        image.image = bpy.data.images["mn_annotations"]

    # node Alpha Over
    alpha_over = mn_compositor.nodes.new("CompositorNodeAlphaOver")
    alpha_over.name = "Alpha Over"
    # Fac
    alpha_over.inputs["Fac"].default_value = 1.0

    # Set locations
    group_output.location = (380.0, 20.0)
    group_input.location = (0.0, 20.0)
    image.location = (0.0, -80.0)
    alpha_over.location = (200.0, -80.0)

    # Set dimensions
    group_output.width, group_output.height = 140.0, 100.0
    group_input.width, group_input.height = 140.0, 100.0
    image.width, image.height = 140.0, 100.0
    alpha_over.width, alpha_over.height = 140.0, 100.0

    # initialize mn_compositor links
    # alpha_over.Image -> group_output.Image
    mn_compositor.links.new(alpha_over.outputs[0], group_output.inputs[0])
    # image.Image -> alpha_over.Image
    mn_compositor.links.new(image.outputs[0], alpha_over.inputs[2])
    # group_input.Image -> alpha_over.Image
    mn_compositor.links.new(group_input.outputs[0], alpha_over.inputs[1])
    return mn_compositor
