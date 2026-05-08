from typing import cast
import bpy
from bpy.types import CompositorNodeTree
from ..nodes.arrange import arrange_tree
from ..nodes.compositor import compositor_node_tree

annotations_image = "mn_annotations"
mn_compositor_node_name = "MN Compositor"


def setup_compositor(scene: bpy.types.Scene):
    node_tree = cast(CompositorNodeTree, scene.compositing_node_group)
    scene.render.use_lock_interface = True
    # add a quick check to see if everything is setup correctly
    if node_tree is not None:
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
        compositor_node_tree()
    # add MN Compositor node to the node tree if not present
    if mn_compositor_node_name not in nodes:
        mn_compositor_node = nodes.new("CompositorNodeGroup")
        mn_compositor_node.node_tree = bpy.data.node_groups[mn_compositor_node_name]
        mn_compositor_node.name = mn_compositor_node_name
    mn_compositor_node = nodes[mn_compositor_node_name]
    # insert MN Compositor right before the Composite node
    # add "Composite" node to node tree if not present
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
