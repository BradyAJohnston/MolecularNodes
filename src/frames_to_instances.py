import bpy

name = "elvator"

collection_frames = bpy.data.collections[name + "_frames"].objects

n_frames = len(collection_frames)

node_group_name = "MOL_" + name + "_frames_to_instances"

new_group = bpy.data.node_groups.get(node_group_name)
if not new_group:
    new_group = bpy.data.node_groups.new(node_group_name, "GeometryNodeTree")

new_group_in = new_group.nodes.new("NodeGroupInput")
new_group_in.location = [-200, 0]

new_group.inputs.new("NodeSocketGeometry", "Atoms")
new_group.inputs.new("NodeSocketCollection", "Frames Collection")
new_group.inputs['Frames Collection'].default_value = collection_frames

new_node = new_group.nodes.new
new_link = new_group.links.new

node_collection_info = new_node("GeometryNodeCollectionInfo")
node_collection_info.location = [0, -100]
node_collection_info.inputs[1].default_value = True
node_collection_info.inputs[2].default_value = True

# initial link of the frames collection from the inputs
new_link(
    new_group_in.outputs['Frames Collection'], 
    node_collection_info.inputs[0]
    )

node_list = list()

for i in range(n_frames):
    node_current = new_node("GeometryNodeGroup")
    node_current.node_tree = bpy.data.node_groups['MOL_utils_split_frame']
    node_current.location = [200, i * -200]
    new_link(
        new_group_in.outputs['Atoms'], 
        node_current.inputs['Atoms']
    )
    
    if i == 0:
        new_link(
            node_collection_info.outputs['Geometry'], 
            node_current.inputs['Frames Remaining']
        )
        node_previous = node_current 
        node_list.append(node_current)
    else:
        new_link(
            node_previous.outputs['Frames Remaining'], 
            node_current.inputs['Frames Remaining']
        )
        node_previous = node_current
        node_list.append(node_current)
    

# create the join geometry node
node_join_geom = new_node("GeometryNodeJoinGeometry")
node_join_geom.location = [500, 0]

node_list.reverse()
for i in range(len(node_list)):
    node = node_list[i]
    new_link(
        node.outputs['Frame Instance'], 
        node_join_geom.inputs[0]
    )


# create the Group Ouptut node
node_group_output = new_node("NodeGroupOutput")
node_group_output.location = [800, 0]

node_group_output.inputs.new("NodeSocketGeometry", "Frame Instances")

new_link(
    node_join_geom.outputs['Geometry'], 
    node_group_output.inputs['Frame Instances']
)
