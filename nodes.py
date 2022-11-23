import bpy
from .tools import property_exists
from .globals import mn_folder
import os

socket_types = {
        'BOOLEAN'  : 'NodeSocketBool', 
        'GEOMETRY' : 'NodeSocketGeometry', 
        'INT'      : 'NodeSocketInt', 
        'MATERIAL' : 'NodeSocketMaterial', 
        'VECTOR'   : 'NodeSocketVector', 
        'STRING'   : 'NodeSocketString', 
        'VALUE'    : 'NodeSocketFloat', 
        'COLLETION': 'NodeSocketCollection', 
        'TEXTURE'  : 'NodeSocketTexture', 
        'COLOR'    : 'NodeSocketColor', 
        'IMAGE'    : 'NodeSocketImage'
    }

def mol_append_node(node_name):
    if bpy.data.node_groups.get(node_name):
        pass
    else:
        before_data = list(bpy.data.node_groups)
        bpy.ops.wm.append(
            directory = os.path.join(
                    mn_folder, 'assets', 'node_append_file.blend' + r'/NodeTree'), 
                    filename = node_name, 
                    link = False
                )   
        new_data = list(filter(lambda d: not d in before_data, list(bpy.data.node_groups)))
    
    return bpy.data.node_groups[node_name]

def mol_base_material():
    """Create MOL_atomic_material. If it already exists, just return the material."""
    mat = bpy.data.materials.get('MOL_atomic_material')
    if not mat:
        mat = bpy.data.materials.new('MOL_atomic_material')
        mat.use_nodes = True
        node_att = mat.node_tree.nodes.new("ShaderNodeAttribute")
        node_att.attribute_name = "Color"
        node_att.location = [-300, 200]
        mat.node_tree.links.new(node_att.outputs['Color'], mat.node_tree.nodes['Principled BSDF'].inputs['Base Color'])
    return mat

def gn_new_group_empty(name = "Geometry Nodes"):
    group = bpy.data.node_groups.get(name)
    # if the group already exists, return it and don't create a new one
    if group:
        return group
    
    # create a new group for this particular name and do some initial setup
    group = bpy.data.node_groups.new(name, 'GeometryNodeTree')
    group.inputs.new('NodeSocketGeometry', "Geometry")
    group.outputs.new('NodeSocketGeometry', "Geometry")
    input_node = group.nodes.new('NodeGroupInput')
    output_node = group.nodes.new('NodeGroupOutput')
    output_node.is_active_output = True
    input_node.select = False
    output_node.select = False
    input_node.location.x = -200 - input_node.width
    output_node.location.x = 200
    group.links.new(output_node.inputs[0], input_node.outputs[0])
    return group

def add_custom_node_group(parent_group, node_name, location = [0,0], width = 200):
    
    mol_append_node(node_name)
    
    node = parent_group.node_group.nodes.new('GeometryNodeGroup')
    node.node_tree = bpy.data.node_groups[node_name]
    
    node.location = location
    node.width = 200 # unsure if width will work TODO check
    
    return node

def add_custom_node_group_to_node(parent_group, node_name, location = [0,0], width = 200):
    
    mol_append_node(node_name)
    
    node = parent_group.nodes.new('GeometryNodeGroup')
    node.node_tree = bpy.data.node_groups[node_name]
    
    node.location = location
    node.width = 200 # unsure if width will work TODO check
    
    return node

def create_starting_node_tree(obj, n_frames = 1, starting_style = "atoms"):
    
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    node_mod = obj.modifiers.get('MolecularNodes')
    if not node_mod:
        node_mod = obj.modifiers.new("MolecularNodes", "NODES")
    obj.modifiers.active = node_mod

    # create a new GN node group, specific to this particular molecule
    node_group = gn_new_group_empty("MOL_" + str(obj.name))
    node_mod.node_group = node_group
    
    # TODO check if can delete this loop
    # ensure the required setup nodes either already exist or append them
    # required_setup_nodes = ['MOL_prop_setup', 'MOL_style_color']
    # if n_frames > 1:
    #     required_setup_nodes = ['MOL_prop_setup', 'MOL_style_color', 'MOL_animate', 'MOL_animate_frames']
    # for node_group in required_setup_nodes:
    #     mol_append_node(node_group)
    
    # move the input and output nodes for the group
    node_input = node_mod.node_group.nodes['Group Input']
    node_input.location = [0, 0]
    node_output = node_mod.node_group.nodes['Group Output']
    node_output.location = [800, 0]
    
    # node_properties = add_custom_node_group(node_group, 'MOL_prop_setup', [0, 0])
    node_colour = add_custom_node_group(node_mod, 'MOL_style_colour', [200, 0])
    
    node_random_colour = node_group.nodes.new("FunctionNodeRandomValue")
    node_random_colour.data_type = 'FLOAT_VECTOR'
    node_random_colour.location = [-60, -200]
    
    node_chain_id = node_group.nodes.new("GeometryNodeInputNamedAttribute")
    node_chain_id.location = [-250, -450]
    node_chain_id.data_type = "INT"
    node_chain_id.inputs['Name'].default_value = "chain_id"
    
    
    
    # create the links between the the nodes that have been established
    link = node_group.links.new
    link(node_input.outputs['Geometry'], node_colour.inputs['Atoms'])
    link(node_colour.outputs['Atoms'], node_output.inputs['Geometry'])
    link(node_random_colour.outputs['Value'], node_colour.inputs['Carbon'])
    link(node_chain_id.outputs[4], node_random_colour.inputs['ID'])
    
    if starting_style == "atoms":
        node_atoms = add_custom_node_group(node_mod, "MOL_style_atoms", location = [500, 0])
        link(node_colour.outputs['Atoms'], node_atoms.inputs['Atoms'])
        link(node_atoms.outputs['Atoms'], node_output.inputs['Geometry'])
        node_atoms.inputs['Material'].default_value = mol_base_material()

    
    # if multiple frames, set up the required nodes for an aniamtion
    if n_frames > 1:
        node_output.location = [1100, 0]
        
        node_animate_frames = add_custom_node_group(node_group, 'MOL_animate_frames', [800, 0])
        
        # node_animate_frames.inputs['Frames Collection'].default_value = col_frames
        node_animate_frames.inputs['Absolute Frame Position'].default_value = True
        
        node_animate = add_custom_node_group(node_group, 'MOL_animate', [550, -300])
        link(node_colour.outputs['Atoms'], node_animate_frames.inputs['Atoms'])
        link(node_animate_frames.outputs['Atoms'], node_output.inputs['Geometry'])
        link(node_animate.outputs['Animate Mapped'], node_animate_frames.inputs[2])


def create_custom_surface(name, n_chains):
    # if a group of this name already exists, just return that instead of making a new one
    group = bpy.data.node_groups.get(name)
    if group:
        return group
    
    # get the node to create a loop from
    looping_node = mol_append_node('MOL_style_surface_single')
    
    
    # create new empty data block
    group = bpy.data.node_groups.new(name, 'GeometryNodeTree')
    
    # loop over the inputs and create an input for each
    for i in looping_node.inputs.values():
        group.inputs.new(socket_types.get(i.type), i.name)
    
    
    
    group.inputs['Selection'].default_value = True
    group.inputs['Selection'].hide_value = True
    group.inputs['Resolution'].default_value = 7
    group.inputs['Radius'].default_value = 1
    group.inputs['Shade Smooth'].default_value = True
    group.inputs['Color By Chain'].default_value = True
    
    # loop over the outputs and create an output for each
    for o in looping_node.outputs.values():
        group.outputs.new(socket_types.get(o.type), o.name)
    
    group.outputs.new(socket_types.get('GEOMETRY'), 'Chain Instances')
    
    # add in the inputs and theo outputs inside of the node
    node_input = group.nodes.new('NodeGroupInput')
    node_input.location = [-300, 0]
    node_output = group.nodes.new('NodeGroupOutput')
    node_output.location = [800, 0]
    
    link = group.links.new
    
    node_input = group.nodes['Group Input']
    # node_input.inputs['Geometry'].name = 'Atoms'
    node_output = group.nodes['Group Output']
    
    node_chain_id = group.nodes.new("GeometryNodeInputNamedAttribute")
    node_chain_id.location = [-250, -450]
    node_chain_id.data_type = "INT"
    node_chain_id.inputs['Name'].default_value = "chain_id"
    
    # for each chain, separate the geometry and choose only that chain, pipe through 
    # a surface node and then join it all back together again
    list_node_surface = []
    height_offset = 300
    for chain in range(n_chains):
        offset = 0 - chain * height_offset
        node_separate = group.nodes.new('GeometryNodeSeparateGeometry')
        node_separate.location = [120, offset]
        
        node_compare = group.nodes.new('FunctionNodeCompare')
        node_compare.data_type = 'INT'
        node_compare.location = [-100, offset]
        node_compare.operation = 'EQUAL'
        
        link(node_chain_id.outputs[4], node_compare.inputs[2])
        node_compare.inputs[3].default_value = chain
        link(node_compare.outputs['Result'], node_separate.inputs['Selection'])
        link(node_input.outputs[0], node_separate.inputs['Geometry'])
        
        node_surface_single = group.nodes.new('GeometryNodeGroup')
        node_surface_single.node_tree = looping_node
        node_surface_single.location = [300, offset]
        
        link(node_separate.outputs['Selection'], node_surface_single.inputs['Atoms'])
        # node_surface_single = add_custom_node_group(group, 'MOL_style_surface_single', [100, offset])
        
        for i in node_surface_single.inputs.values():
            if i.type != 'GEOMETRY':
                link(node_input.outputs[i.name], i)
        
        
        list_node_surface.append(node_surface_single)
    
    # create join geometry, and link the nodes in reverse order
    node_join_geometry = group.nodes.new('GeometryNodeJoinGeometry')
    node_join_geometry.location = [500, 0]
    
    node_instance = group.nodes.new('GeometryNodeGeometryToInstance')
    node_instance.location = [500, -600]
    
    node_join_volume = group.nodes.new('GeometryNodeJoinGeometry')
    node_join_volume.location = [500, -300]
    
    list_node_surface.reverse()
    
    # TODO: turn these into instances instead of just joining geometry
    for n in list_node_surface:
        link(n.outputs['Surface'], node_join_geometry.inputs['Geometry'])
        link(n.outputs['Surface'], node_instance.inputs['Geometry'])
        link(n.outputs['Volume'], node_join_volume.inputs['Geometry'])
    
    
    
    link(node_join_geometry.outputs['Geometry'], node_output.inputs[0])
    link(node_instance.outputs['Instances'], node_output.inputs['Chain Instances'])
    link(node_join_volume.outputs['Geometry'], node_output.inputs['Volume'])
    
    return group

def rotation_matrix(node_group, mat_rot, mat_trans, location = [0,0]):
    
    node_utils_rot = mol_append_node('MOL_utils_rotation_matrix')
    
    node = node_group.nodes.new('GeometryNodeGroup')
    node.node_tree = node_utils_rot
    node.location = location
    
    for rot in range(3):
        for value in range(3):
            node.inputs[rot].default_value[value] = mat_rot[rot, value]
    
    for value in range(3):
        node.inputs[3].default_value[value] = mat_trans[value]
    
    return node