import bpy

# all that should be passed into the actual function
# these aren't defined in the script but the variables are set up by 
# serpens when the addon is compiled
# 
# # molecule_name = "elvator"
# style = "style_atoms"

node_group_name = "MOL_" + molecule_name + "_realise_frames_" + style

# try and get the frames collection, throw an error if it doesn't exist
collection_frames = bpy.data.collections.get(molecule_name + "_frames")
assert collection_frames, "Collection " + molecule_name + "_frames does not exist."

n_frames = len(collection_frames.objects)

def node_iterate_join(n_frames, split_node_group_name = 'MOL_utils_split_frames', node_height = 200):
    
    iterated_node_group = bpy.data.node_groups[split_node_group_name]    
    
    # list of newly created nodes, that can then be iterated over and 
    # linked to the join geometry node after creation
    node_list = list()

    for i in range(n_frames):
        node_current = new_node("GeometryNodeGroup")
        node_current.node_tree = iterated_node_group
        node_current.location = [200, i * -node_height]
        new_link(
            node_group_in.outputs['Atoms'], 
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
    
    return node_list

# define the properties of the different nodes that will require iterating
# can probably automate this part in the future, but for now just going to do this
dict_props_ribbon = {
    "Smooth Curve"      :    {"default": 0, "name": "Int"}, 
    "Radius Ribbon"     :    {"default": 3, "name": "Float"}, 
    "Resolution Ribbon" :    {"default": 8, "name": "Int"}, 
    "Profile Ribbon"    :                  {"name": "Geometry"}, 
    "Shade Smooth"      : {"default": True, "name": "Bool"}, 
    "Material"          :                  {"name": "Material"}
}

dict_props_surface  = {
    "Ignore Hydrogens" : {"default" :  True, "name" : "Bool"},
    "Radius"           : {"default" :   1.5, "name" : "Float"},
    "Voxel Size"       : {"default" :   0.2, "name" : "Float"},
    "Threshold"        : {"default" :   0.1, "name" : "Float"},
    "Adaptivity"       : {"default" :     0, "name" : "Float"},
    "Shade Smooth"     : {"default" :  True, "name" : "Bool"},
    "Merge Distance"   : {"default" : 0.001, "name" : "Float"},
    "Smoothing"        : {"default" :     0, "name" : "Int"},
    "Color by Chain"   : {"default" :  True, "name" : "Bool"},
    "Material"         :                    {"name" : "Material"}
}

dict_props_atoms_EEVEE = {
    "Subdivisions"  : {"defualt" : 1,         "name" : "Int"},
    "Scale Radii"   : {"defualt" : 1,         "name" : "Float"},
    "Custom Atom"   :                        {"name" : "Geometry"},
    "Atom Rotation" : {"defualt" : [0, 0, 0], "name" : "Vector"},
    "Shade Smooth"  : {"defualt" : True,      "name" : "Bool"},
    "Material"      :                        {"name" : "Material"}
}

dict_props_atoms = {}

dict_style = {
    "style_atoms" :       {"node" : "MOL_utils_split_frame", "dict" : dict_props_atoms}, 
    "style_atoms_EEVEE" : {"node" : "MOL_utils_split_frame_atoms_EEVEE", "dict" : dict_props_atoms_EEVEE},
    "style_ribbon" :      {"node" : "MOL_utils_split_frame_ribbon", "dict" : dict_props_ribbon},
    "style_surface" :     {"node" : "MOL_utils_split_frame_surface", "dict" : dict_props_surface}
}


# Attempt to find the node group
new_group = bpy.data.node_groups.get(node_group_name)
# If the group doesn't exist, create it and set up the required nodes
if not new_group:
    new_group = bpy.data.node_groups.new(node_group_name, "GeometryNodeTree")

    node_group_in = new_group.nodes.new("NodeGroupInput")
    node_group_in.location = [-200, 0]

    new_group.inputs.new("NodeSocketGeometry", "Atoms")
    new_group.inputs.new("NodeSocketCollection", "Frames Collection")

    
    # add boolean for in the input to switch between computing single and all frames
    new_group.inputs.new("NodeSocketBool", "Single Frame")
    new_group.inputs['Single Frame'].default_value = True
    # new_group.inputs['Frames Collection'].default_value = collection_frames.name

    new_node = new_group.nodes.new
    new_link = new_group.links.new

    node_collection_info = new_node("GeometryNodeCollectionInfo")
    node_collection_info.location = [0, -100]
    node_collection_info.inputs[1].default_value = True
    node_collection_info.inputs[2].default_value = True

    # initial link of the frames collection from the inputs
    new_link(
        node_group_in.outputs['Frames Collection'], 
        node_collection_info.inputs[0]
        )

    
    # get the relevant dictionary for the given style that will be used
    # inside of the following loops
    loop_dict = dict_style.get(style).get("dict")
    
    # iterate through and create a copy of the specified node for each frame
    # doe some initial linking to each other in the chain, however the style-specific
    # linking needs to be done in a second loop
    node_list = node_iterate_join(
        n_frames = n_frames,  
        split_node_group_name = dict_style.get(style).get("node"), # get the node name from the dictionary to iterate
        node_height = 200
    )


    # Setup the properties for the given style type
    for prop in loop_dict.keys():
        new_group.inputs.new("NodeSocket" + loop_dict[prop].get("name"), prop)
        try:
            new_group.inputs[prop].default_value = loop_dict.get(prop).get("default")
        except:
            next

    for node in node_list:
        for prop in loop_dict.keys():
            new_link(
                node_group_in.outputs[prop], 
                node.inputs[prop]
            )


    # create the join geometry node
    node_join_geom = new_node("GeometryNodeJoinGeometry")
    node_join_geom.location = [500, 0]

    # When linking nodes to a "Join Geometry" node, it auto-joins them in 
    # the reverse wrong oder. First we reverse the order of the list, and then 
    # go through the list and join all of the nodes to the "Join Geometry" node 
    node_list.reverse()
    for i in range(len(node_list)):
        node = node_list[i]
        new_link(
            node.outputs['Frame Instance'], 
            node_join_geom.inputs[0]
        )


    # create the Group Ouptut node
    node_group_out = new_node("NodeGroupOutput")
    node_group_out.location = [1000, 0]
    new_group.outputs.new("NodeSocketGeometry", "Frame Instances")

    
    # create the switch for changing between single frame
    # and all frames
    node_switch = new_node("GeometryNodeSwitch")
    node_switch.location = [750, 0]

    # link the boolean to the geometry switch
    new_link(
        node_group_in.outputs['Single Frame'], 
        node_switch.inputs[1]
    )

    # link the "Join Geometry" node to the swithc
    new_link(
        node_join_geom.outputs[0], 
        node_switch.inputs[14]
    )
    # re-reverse the list to easily access the first node created
    node_list.reverse()
    first_node = node_list[0]

    # link the first 
    new_link(
        first_node.outputs[0], 
        node_switch.inputs[15]
    )
    

    # link the switch to the output of the node group
    new_link(
        node_switch.outputs[6], 
        node_group_out.inputs[0]
    )
