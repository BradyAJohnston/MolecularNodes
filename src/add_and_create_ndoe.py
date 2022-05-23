import bpy

# Get the currently active GN tree
obj = bpy.context.active_object
node_tree = obj.modifiers.active


# new_node = new_node
# Add the new node group, will have been appended to the scene previously by the 'load_node' function
bpy.ops.node.add_node(
    type = "GeometryNodeGroup", 
    use_transform = True, 
    settings = [
        {"name":"node_tree", "value":"bpy.data.node_groups['" + new_node + "']"}
        ]
        )

# select the just-added node, and make it currently being transformed by the user
node = bpy.context.selected_nodes[0]
bpy.ops.transform.translate('INVOKE_DEFAULT')

