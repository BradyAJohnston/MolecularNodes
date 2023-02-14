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
    """Append MOL_atomic_material to the .blend file it it doesn't already exist, and return that material."""
    mat = bpy.data.materials.get('MOL_atomic_material')
    
    if not mat:
        mat = bpy.ops.wm.append(
            directory=os.path.join(
                mn_folder, 'assets', 'node_append_file.blend' + r'/Material'
            ), 
            filename='MOL_atomic_material', 
            link=False
        )
    
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

def create_starting_node_tree(obj, coll_frames, starting_style = "atoms"):
    
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
    node_input = node_mod.node_group.nodes[bpy.app.translations.pgettext_data("Group Input",)]
    node_input.location = [0, 0]
    node_output = node_mod.node_group.nodes[bpy.app.translations.pgettext_data("Group Output",)]
    node_output.location = [800, 0]
    
    # node_properties = add_custom_node_group(node_group, 'MOL_prop_setup', [0, 0])
    node_colour = add_custom_node_group(node_mod, 'MOL_style_color', [200, 0])
    
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
    
    styles = ['MOL_style_atoms_cycles', 'MOL_style_ribbon_protein', 'MOL_style_ball_and_stick']
    
    # if starting_style == "atoms":
    
    node_style = add_custom_node_group(node_mod, styles[starting_style], location = [500, 0])
    link(node_colour.outputs['Atoms'], node_style.inputs['Atoms'])
    link(node_style.outputs[0], node_output.inputs['Geometry'])
    node_style.inputs['Material'].default_value = mol_base_material()

    
    # if multiple frames, set up the required nodes for an aniamtion
    if coll_frames:
        node_output.location = [1100, 0]
        node_style.location = [800, 0]
        
        node_animate_frames = add_custom_node_group_to_node(node_group, 'MOL_animate_frames', [500, 0])
        node_animate_frames.inputs['Frames'].default_value = coll_frames
        
        # node_animate_frames.inputs['Absolute Frame Position'].default_value = True
        
        node_animate = add_custom_node_group_to_node(node_group, 'MOL_animate_value', [500, -300])
        link(node_colour.outputs['Atoms'], node_animate_frames.inputs['Atoms'])
        link(node_animate_frames.outputs['Atoms'], node_style.inputs['Atoms'])
        link(node_animate.outputs['Animate 0..1'], node_animate_frames.inputs['Animate 0..1'])


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
    
    node_input = group.nodes[bpy.app.translations.pgettext_data("Group Input",)]
    # node_input.inputs['Geometry'].name = 'Atoms'
    node_output = group.nodes[bpy.app.translations.pgettext_data("Group Output",)]
    
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

def rotation_matrix(node_group, mat, location = [0,0], world_scale = 0.01):
    """Add a Rotation & Translation node from a 3x4 matrix.

    Args:
        node_group (_type_): Parent node group to add this new node to.
        mat (_type_): 3x4 rotation & translation matrix
        location (list, optional): Position to add the node in the node tree. Defaults to [0,0].
        world_scale(float, optional): Scaling factor for the world. Defaults to 0.01.
    Returns:
        _type_: Newly created node tree.
    """
    from scipy.spatial.transform import Rotation as R
    
    node_utils_rot = mol_append_node('MOL_utils_rot_trans')
    
    node = node_group.nodes.new('GeometryNodeGroup')
    node.node_tree = node_utils_rot
    node.location = location
    
    # calculate the euler rotation from the rotation matrix
    rotation = R.from_matrix(mat[:3, :3]).as_euler('xyz')
    
    # set the values for the node that was just created
    # set the euler rotation values
    for i in range(3):
        node.inputs[0].default_value[i] = rotation[i]
    # set the translation values
    for i in range(3):
        node.inputs[1].default_value[i] = mat[:3, 3:][i] * world_scale
        
    return node

def chain_selection(node_name, input_list, attribute, starting_value = 0, label_prefix = ""):
    """
    Given a an input_list, will create a node which takes an Integer input, 
    and has a boolean tick box for each item in the input list. The outputs will
    be the resulting selection and the inversion of the selection.
    Can contain a prefix for the resulting labels. Mostly used for constructing 
    chain selections when required for specific proteins.
    """
    # just reutn the group name early if it already exists
    group = bpy.data.node_groups.get(node_name)
    if group:
        return group
    
    # get the active object, might need to change to taking an object as an input
    # and making it active isntead, to be more readily applied to multiple objects
    obj = bpy.context.active_object
    # try to get the Molecular Nodes modifier and select it, if not create one and select it
    node_mod = obj.modifiers.get('MolecularNodes')
    if not node_mod:
        node_mod = obj.modifiers.new("MolecularNodes", "NODES")
    
    obj.modifiers.active = node_mod
    
    
    # link shortcut for creating links between nodes
    link = node_mod.node_group.links.new
    # create the custom node group data block, where everything will go
    # also create the required group node input and position it
    chain_group = bpy.data.node_groups.new(node_name, "GeometryNodeTree")
    chain_group_in = chain_group.nodes.new("NodeGroupInput")
    chain_group_in.location = [-200, 0]
    # create a named attribute node that gets the chain_number attribute
    # and use this for the selection algebra that happens later on
    chain_number_node = chain_group.nodes.new("GeometryNodeInputNamedAttribute")
    chain_number_node.data_type = 'INT'
    chain_number_node.location = [-200, 200]
    chain_number_node.inputs[0].default_value = attribute
    chain_number_node.outputs.get('Attribute')
    # create a boolean input for the group for each item in the list
    for chain_name in input_list: 
        # create a boolean input for the name, and name it whatever the chain chain name is
        chain_group.inputs.new("NodeSocketBool", str(label_prefix) + str(chain_name))
    # shortcut for creating new nodes
    new_node = chain_group.nodes.new
    # distance horizontally to space all of the created nodes
    node_sep_dis = 180
    counter = 0
    for chain_name in input_list:
        current_node = chain_group.nodes.new("GeometryNodeGroup")
        current_node.node_tree = mol_append_node('MOL_utils_bool_chain')
        current_node.location = [counter * node_sep_dis, 200]
        current_node.inputs["number_matched"].default_value = counter + starting_value
        group_link = chain_group.links.new
        # link from the the named attribute node chain_number into the other inputs
        if counter == 0:
            # for some reason, you can't link with the first output of the named attribute node. Might
            # be a bug, which might be changed later, so I am just going through a range of numbers for 
            # the named attribute node outputs, to link whatever it ends up being. Dodgy I know. 
            # TODO revisit this and see if it is fixed and clean up code
            for i in range(5):
                try:
                    group_link(chain_number_node.outputs[i], current_node.inputs['number_chain_in'])
                except:
                    pass
        group_link(chain_group_in.outputs[counter], current_node.inputs["bool_include"])
        if counter > 0:
            group_link(previous_node.outputs['number_chain_out'], current_node.inputs['number_chain_in'])
            group_link(previous_node.outputs['bool_chain_out'], current_node.inputs['bool_chain_in'])
        previous_node = current_node
        counter += 1
    chain_group_out = chain_group.nodes.new("NodeGroupOutput")
    chain_group_out.location = [(counter + 1) * node_sep_dis, 200]
    chain_group.outputs.new("NodeSocketBool", "Selection")
    chain_group.outputs.new("NodeSocketBool", "Inverted")
    group_link(current_node.outputs['bool_chain_out'], chain_group_out.inputs['Selection'])
    bool_math = chain_group.nodes.new("FunctionNodeBooleanMath")
    bool_math.location = [counter * node_sep_dis, 50]
    bool_math.operation = "NOT"
    group_link(current_node.outputs['bool_chain_out'], bool_math.inputs[0])
    group_link(bool_math.outputs[0], chain_group_out.inputs['Inverted'])
    # create an empty node group group inside of the node tree
    # link the just-created custom node group data to the node group in the tree
    # new_node_group = node_mod.node_group.nodes.new("GeometryNodeGroup")
    # new_node_group.node_tree = bpy.data.node_groups[chain_group.name]
    # resize the newly created node to be a bit wider
    # node_mod.node_group.nodes[-1].width = 200
    # the chain_id_list and output_name are passed in from the operator when it is called
    # these are custom properties that are associated with the object when it is initial created
    return chain_group

def chain_color(node_name, input_list, label_prefix = "Chain "):
    """
    Given the input list of chain names, will create a node group which uses
    the chain_id named attribute to manually set the colours for each of the chains.
    """
    
    import random
    
    # get the active object, might need to change to taking an object as an input
    # and making it active isntead, to be more readily applied to multiple objects
    obj = bpy.context.active_object
    # try to get the Molecular Nodes modifier and select it, if not create one and select it
    node_mod = obj.modifiers.get('MolecularNodes')
    if not node_mod:
        node_mod = obj.modifiers.new("MolecularNodes", "NODES")
    
    obj.modifiers.active = node_mod
    
    
    # create the custom node group data block, where everything will go
    # also create the required group node input and position it
    chain_group = bpy.data.node_groups.new(node_name, "GeometryNodeTree")
    node_input = chain_group.nodes.new("NodeGroupInput")
    node_input.location = [-200, 0]
    
    
    # link shortcut for creating links between nodes
    link = chain_group.links.new
    
    
    # create a named attribute node that gets the chain_number attribute
    # and use this for the selection algebra that happens later on
    chain_number_node = chain_group.nodes.new("GeometryNodeInputNamedAttribute")
    chain_number_node.data_type = 'INT'
    chain_number_node.location = [-200, 400]
    chain_number_node.inputs[0].default_value = 'chain_id'
    chain_number_node.outputs.get('Attribute')
    
    # shortcut for creating new nodes
    new_node = chain_group.nodes.new
    # distance horizontally to space all of the created nodes
    node_sep_dis = 180
    counter = 0
    
    for chain_name in input_list:
        offset = counter * node_sep_dis
        current_chain = str(label_prefix) + str(chain_name)
         # node compare inputs 2 & 3
        node_compare = new_node('FunctionNodeCompare')
        node_compare.data_type = 'INT'
        node_compare.location = [offset, 100]
        node_compare.operation = 'EQUAL'
        
        node_compare.inputs[3].default_value = counter
        
        # link the named attribute to the compare
        link(chain_number_node.outputs[4], node_compare.inputs[2])
        
        node_color = new_node('GeometryNodeSwitch')
        node_color.input_type = 'RGBA'
        node_color.location = [offset, -100]
        
        # create an input for this chain
        chain_group.inputs.new("NodeSocketColor", current_chain)
        chain_group.inputs[current_chain].default_value = [random.random(), random.random(), random.random(), 1]
        # switch input colours 10 and 11
        link(node_input.outputs[current_chain], node_color.inputs[11])
        link(node_compare.outputs['Result'], node_color.inputs['Switch'])
        
        
        if counter > 0:
            link(node_color_previous.outputs[4], node_color.inputs[10])
        
        node_color_previous = node_color
        counter += 1
    chain_group.outputs.new("NodeSocketColor", "Color")
    node_output = chain_group.nodes.new("NodeGroupOutput")
    node_output.location = [offset, 200]
    link(node_color.outputs[4], node_output.inputs['Color'])
    
    return chain_group

def resid_multiple_selection(node_name, input_resid_string):
    """
    Returns a node group that takes an integer input and creates a boolean 
    tick box for each item in the input list. Outputs are the selected 
    residues and the inverse selection. Used for constructing chain 
    selections in specific proteins.
    """
        
    #print(f'recieved input: {input_resid_string}')
    # do a cleanning of input string to allow fuzzy input from users
    for c in ";/+ .":
        if c in input_resid_string:
            input_resid_string=input_resid_string.replace(c, ',')

    for c in "_=:":
        if c in input_resid_string:
            input_resid_string=input_resid_string.replace(c, '-')

    #print(f'fixed input:{input_resid_string}')

    # parse input_resid_string into sub selecting string list
    sub_list=[item for item in input_resid_string.split(',') if item]
    
    # distance vertical to space all of the created nodes
    node_sep_dis = -100

    # get the active object, might need to change to taking an object as an input
    # and making it active isntead, to be more readily applied to multiple objects
    obj = bpy.context.active_object
    # try to get the Molecular Nodes modifier and select it, if not create one and select it
    node_mod = obj.modifiers.get('MolecularNodes')
    if not node_mod:
        node_mod = obj.modifiers.new("MolecularNodes", "NODES")
    
    obj.modifiers.active = node_mod
    
    # create the custom node group data block, where everything will go
    # also create the required group node input and position it
    residue_id_group = bpy.data.node_groups.new(node_name, "GeometryNodeTree")
    residue_id_group_in = residue_id_group.nodes.new("NodeGroupInput")
    residue_id_group_in.location = [0, node_sep_dis * len(sub_list)/2]
    
    group_link = residue_id_group.links.new
    new_node = residue_id_group.nodes.new
    
    for residue_id_index,residue_id in enumerate(sub_list):
        
        if '-' in residue_id:
            [resid_start, resid_end] = residue_id.split('-')[:2]
            # set two new inputs
            residue_id_group.inputs.new("NodeSocketInt",'res_id_start').default_value = int(resid_start)
            residue_id_group.inputs.new("NodeSocketInt",'res_id_end').default_value = int(resid_end)
        else:
            # set a new input and set the resid
            residue_id_group.inputs.new("NodeSocketInt",'res_id').default_value = int(residue_id)
        
    # set a counter for MOL_sel_res_id* nodes
    counter=0
    for residue_id_index,residue_id in enumerate(sub_list):
        
        # add an new node of MOL_sel_res_id or MOL_sek_res_id_range
        current_node = new_node("GeometryNodeGroup")

        # add an bool_math block 
        bool_math = new_node("FunctionNodeBooleanMath")
        bool_math.location = [400,(residue_id_index+1) * node_sep_dis]
        bool_math.operation = "OR"

        if '-' in residue_id:
            # a residue range
            current_node.node_tree = mol_append_node('MOL_sel_res_id_range')
            
            group_link(residue_id_group_in.outputs[counter], current_node.inputs[0])
            counter+=1
            
            group_link(residue_id_group_in.outputs[counter], current_node.inputs[1])
            
        else:
            # create a node
            current_node.node_tree = mol_append_node('MOL_sel_res_id')
            # link the input of MOL_sel_res_id
            #print(f'counter={counter} of {residue_id}')
            group_link(residue_id_group_in.outputs[counter], current_node.inputs[0])
        
        counter+=1
        
        # set the coordinates
        current_node.location = [200,(residue_id_index+1) * node_sep_dis]
    
        if residue_id_index == 0:
            # link the first residue selection to the first input of its OR block
            group_link(current_node.outputs['Selection'],bool_math.inputs[0])
        else:
            # if it is not the first residue selection, link the output to the previous or block
            group_link(current_node.outputs['Selection'], previous_bool_node.inputs[1])
        
            # link the ouput of previous OR block to the current OR block
            group_link(previous_bool_node.outputs[0], bool_math.inputs[0])
        previous_bool_node = bool_math

    # add a output block
    residue_id_group_out = new_node("NodeGroupOutput")
    residue_id_group_out.location = [800,(residue_id_index + 1) / 2 * node_sep_dis]
    residue_id_group.outputs.new("NodeSocketBool", "Selection")
    residue_id_group.outputs.new("NodeSocketBool", "Inverted")
    group_link(previous_bool_node.outputs[0], residue_id_group_out.inputs['Selection'])
    invert_bool_math = new_node("FunctionNodeBooleanMath")
    invert_bool_math.location = [600,(residue_id_index+1)/ 3 * 2 * node_sep_dis]
    invert_bool_math.operation = "NOT"
    group_link(previous_bool_node.outputs[0], invert_bool_math.inputs[0])
    group_link(invert_bool_math.outputs[0], residue_id_group_out.inputs['Inverted'])
    return residue_id_group