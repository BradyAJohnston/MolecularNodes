import bpy
from ..nodes.arrange import arrange_tree
from ..nodes.compositor import default_5x_compositor_node_tree, mn_compositor_node_tree

annotations_image = "mn_annotations"
mn_compositor_node_name = "MN Compositor"


def setup_compositor(scene: bpy.types.Scene):
    if bpy.app.version >= (5, 0, 0):
        node_tree = scene.compositing_node_group
        use_nodes = True
    else:
        node_tree = scene.node_tree
        use_nodes = scene.use_nodes
    # lock interface when rendering
    scene.render.use_lock_interface = True
    # add a quick check to see if everything is setup correctly
    if node_tree:
        if (
            use_nodes
            and mn_compositor_node_name in node_tree.nodes
            and node_tree.nodes[mn_compositor_node_name].inputs["Image"].is_linked
            and node_tree.nodes[mn_compositor_node_name].outputs["Image"].is_linked
        ):
            # MN Compositor node is present and both inputs and output are linked
            return
    else:
        # no node tree, create one
        if bpy.app.version >= (5, 0, 0):
            # Staring 5.x compositor node trees can be re-used
            # Technically we don't have to create a new one if we import from
            # assets, but this is the safest when generating from code
            scene.compositing_node_group = default_5x_compositor_node_tree()
            node_tree = scene.compositing_node_group
        else:
            scene.use_nodes = True
            node_tree = scene.node_tree
    # setup the compositor node tree
    nodes = node_tree.nodes
    links = node_tree.links
    # create placeholder annotation image if not present
    if annotations_image not in bpy.data.images:
        bpy.data.images.new(annotations_image, 1, 1)
    # add MN Compositor node group to data block if not present
    if mn_compositor_node_name not in bpy.data.node_groups:
        mn_compositor_node_tree()
    # add MN Compositor node to the node tree if not present
    if mn_compositor_node_name not in nodes:
        mn_compositor_node = nodes.new("CompositorNodeGroup")
        mn_compositor_node.node_tree = bpy.data.node_groups[mn_compositor_node_name]
        mn_compositor_node.name = mn_compositor_node_name
    mn_compositor_node = nodes[mn_compositor_node_name]
    # insert MN Compositor right before the Composite node
    # add "Composite" node to node tree if not present
    if bpy.app.version >= (5, 0, 0):
        if "Group Output" not in nodes:
            nodes.new(type="NodeGroupOutput")
        output_node = nodes["Group Output"]
    else:
        if "Composite" not in nodes:
            nodes.new(type="CompositorNodeComposite")
        output_node = nodes["Composite"]
    # add "Render Layers" node to node tree if not present
    if "Render Layers" not in nodes:
        nodes.new(type="CompositorNodeRLayers")
    render_layers_node = nodes["Render Layers"]
    if output_node.inputs["Image"].is_linked:
        # Composite node input is linked - ensure it is from MN Compositor
        from_node = output_node.inputs["Image"].links[0].from_node
        if from_node.name != mn_compositor_node_name:
            # Composite node input is not MN Compositor node, re-link correctly
            from_socket = output_node.inputs["Image"].links[0].from_socket
            links.new(from_socket, mn_compositor_node.inputs["Image"])
            links.new(mn_compositor_node.outputs["Image"], output_node.inputs["Image"])
    else:
        # Composite node input is not linked
        from_socket = render_layers_node.outputs["Image"]
        links.new(from_socket, mn_compositor_node.inputs["Image"])
        links.new(mn_compositor_node.outputs["Image"], output_node.inputs["Image"])
    # link MN Compositor node input if not linked
    if not mn_compositor_node.inputs["Image"].is_linked:
        from_socket = render_layers_node.outputs["Image"]
        links.new(from_socket, mn_compositor_node.inputs["Image"])
    # link to Viewer node if present
    if "Viewer" in nodes:
        viewer = nodes["Viewer"]
        links.new(mn_compositor_node.outputs["Image"], viewer.inputs["Image"])
    # arrange the composite node tree
    arrange_tree(node_tree, add_group_input=False)
