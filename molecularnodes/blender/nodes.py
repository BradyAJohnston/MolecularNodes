import bpy
import os
import numpy as np
import math
import warnings
from .. import assembly
from .. import color
from .. import pkg
from ..blender import obj

socket_types = {
        'BOOLEAN'   : 'NodeSocketBool', 
        'GEOMETRY'  : 'NodeSocketGeometry', 
        'INT'       : 'NodeSocketInt', 
        'MATERIAL'  : 'NodeSocketMaterial', 
        'VECTOR'    : 'NodeSocketVector', 
        'STRING'    : 'NodeSocketString', 
        'VALUE'     : 'NodeSocketFloat', 
        'COLLECTION': 'NodeSocketCollection', 
        'TEXTURE'   : 'NodeSocketTexture', 
        'COLOR'     : 'NodeSocketColor', 
        'IMAGE'     : 'NodeSocketImage'
    }

# current implemented representations
styles_mapping = {
    "presets"       : "MN_style_presets",
    'preset_1'      : ".MN_style_preset_1",
    'preset_2'      : ".MN_style_preset_2",
    'preset_3'      : ".MN_style_preset_3",
    'preset_4'      : ".MN_style_preset_4",
    'atoms'         : 'MN_style_spheres',
    'spheres'       : 'MN_style_spheres',
    'vdw'           : 'MN_style_spheres',
    'sphere'        : 'MN_style_spheres',
    'cartoon'       : 'MN_style_cartoon',
    'ribbon'        : 'MN_style_ribbon',
    'surface'       : 'MN_style_surface',
    'ball_and_stick': 'MN_style_ball_and_stick',
    'ball+stick'    : 'MN_style_ball_and_stick',
    'oxdna'         : 'MN_oxdna_style_ribbon', 
    "density_surface": "MN_density_style_surface", 
    "density_wire" : "MN_density_style_wire"
}

STYLE_ITEMS = (
        ('presets', 'Presets', 'A pre-made combination of different styles'),
        ("spheres", "Spheres", "Space-filling atoms style."), 
        ("surface", "Surface", "Solvent-accsible surface."),
        ("cartoon", "Cartoon", "Secondary structure cartoons"), 
        ("ribbon", "Ribbon", "Continuous backbone ribbon."), 
        ("ball_and_stick", "Ball and Stick", "Spheres for atoms, sticks for bonds")
    )

bpy.types.Scene.MN_import_style = bpy.props.EnumProperty(
    name = "Style", 
    description = "Default style for importing molecules.", 
    items = STYLE_ITEMS
)


MN_DATA_FILE = os.path.join(pkg.ADDON_DIR, 'assets', 'MN_data_file.blend')

def inputs(node):
    items = {}
    for item in node.interface.items_tree:
        if item.item_type == "SOCKET":
            if item.in_out == "INPUT":
                items[item.name] = item 
    return items

def outputs(node):
    items = {}
    for item in node.interface.items_tree:
        if item.item_type == "SOCKET":
            if item.in_out == "OUTPUT":
                items[item.name] = item
    return items

def get_output_type(node, type = "INT"):
    for output in node.outputs:
        if output.type == type:
            return output

def set_selection(group, node, selection):
    pos = node.location
    pos = [pos[0] - 200, pos[1] - 200]
    selection.location = pos
    group.links.new(selection.outputs[0], node.inputs['Selection'])
    
    return selection

def add_selection(group, sel_name, input_list, attribute = 'chain_id'):
    style = style_node(group)
    sel_node = add_custom(group, chain_selection('selection', input_list, attribute=attribute).name)
    
    set_selection(group, style, sel_node)
    return sel_node

def get_output(group):
    return group.nodes[bpy.app.translations.pgettext_data("Group Output",)]

def get_input(group):
    return group.nodes[bpy.app.translations.pgettext_data("Group Input",)]

def get_mod(object):
    node_mod = object.modifiers.get('MolecularNodes')
    if not node_mod:
        node_mod = object.modifiers.new("MolecularNodes", "NODES")
    object.modifiers.active = node_mod
    return node_mod

def format_node_name(name):
    "Formats a node's name for nicer printing."
    return name.strip("MN_").replace("_", " ").title()

def  get_nodes_last_output(group):
    output = get_output(group)
    last = output.inputs[0].links[0].from_node
    return last, output

def previous_node(node):
    "Get the node which is the first connection to the first input of this node"
    prev = node.inputs[0].links[0].from_node
    return prev

def style_node(group):
    prev = previous_node(get_output(group))
    is_style_node = ("style" in prev.name)
    while not is_style_node:
        print(prev.name)
        prev = previous_node(prev)
        is_style_node = ("style" in prev.name)
    return prev

def get_style_node(object):
    "Walk back through the primary node connections until you find the first style node"
    group = object.modifiers['MolecularNodes'].node_group
    return style_node(group)

def star_node(group):
    prev = previous_node(get_output(group))
    is_star_node = ("MN_starfile_instances" in prev.name)
    while not is_star_node:
        print(prev.name)
        prev = previous_node(prev)
        is_star_node = ("MN_starfile_instances" in prev.name)
    return prev

def get_star_node(object):
    "Walk back through the primary node connections until you find the first style node"
    group = object.modifiers['MolecularNodes'].node_group
    return star_node(group)

def get_color_node(object):
    "Walk back through the primary node connections until you find the first style node"
    group = object.modifiers['MolecularNodes'].node_group
    for node in group.nodes:
        if node.name == "MN_color_attribute_random":
            return node

def insert_last_node(group, node):
    last, output = get_nodes_last_output(group)
    link = group.links.new
    location = output.location
    output.location = [location[0] + 300, location[1]]
    node.location = [location[0] - 300, location[1]]
    link(last.outputs[0], node.inputs[0])
    link(node.outputs[0], output.inputs[0])

def realize_instances(obj):
    group = obj.modifiers['MolecularNodes'].node_group
    realize = group.nodes.new('GeometryNodeRealizeInstances')
    insert_last_node(group, realize)

def append(node_name, link = False):
    node = bpy.data.node_groups.get(node_name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if not node or link:
            bpy.ops.wm.append(
                'EXEC_DEFAULT',
                directory = os.path.join(MN_DATA_FILE, 'NodeTree'), 
                filename = node_name, 
                link = link
            )
    
    return bpy.data.node_groups[node_name]

def MN_base_material():
    """
    Append MN_atomic_material to the .blend file it it doesn't already exist, 
    and return that material.
    """
    
    mat_name = 'MN_atomic_material'
    mat = bpy.data.materials.get(mat_name)
    
    if not mat:
        print('appending material')
        bpy.ops.wm.append(
            directory = os.path.join(MN_DATA_FILE, 'Material'), 
            filename = 'MN_atomic_material', 
            link = False
        )
    
    return bpy.data.materials[mat_name]

def new_group(name = "Geometry Nodes", geometry = True, fallback=True):
    group = bpy.data.node_groups.get(name)
    # if the group already exists, return it and don't create a new one
    if group and fallback:
        return group
    
    # create a new group for this particular name and do some initial setup
    group = bpy.data.node_groups.new(name, 'GeometryNodeTree')
    input_node = group.nodes.new('NodeGroupInput')
    output_node = group.nodes.new('NodeGroupOutput')
    input_node.location.x = -200 - input_node.width
    output_node.location.x = 200
    if geometry:
        group.interface.new_socket('Geometry', in_out='INPUT', socket_type='NodeSocketGeometry')
        group.interface.new_socket('Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
        group.links.new(output_node.inputs[0], input_node.outputs[0])
    return group

def assign_material(node, material = 'default'):
    material_socket = node.inputs.get('Material')
    if material_socket:
        if not material:
            pass
        elif material == "default":
            material_socket.default_value = MN_base_material()
        else:
            material_socket.default_value = material

def add_node(node_name, label: str = '', show_options = False, material = "default"):
    # intended to be called upon button press in the node tree
    
    prev_context = bpy.context.area.type
    bpy.context.area.type = 'NODE_EDITOR'
    # actually invoke the operator to add a node to the current node tree
    # use_transform=True ensures it appears where the user's mouse is and is currently 
    # being moved so the user can place it where they wish
    bpy.ops.node.add_node(
        'INVOKE_DEFAULT', 
        type='GeometryNodeGroup', 
        use_transform=True
        )
    bpy.context.area.type = prev_context
    node = bpy.context.active_node
    node.node_tree = bpy.data.node_groups[node_name]
    node.width = 200.0
    node.show_options = show_options
    
    if label == '':
        node.label = format_node_name(node_name)
    else:
        node.label = label
    node.name = node_name
    
    # if added node has a 'Material' input, set it to the default MN material
    assign_material(node, material=material)

def add_custom(
        group, 
        name, 
        location = [0,0], 
        width = 200,
        material = "default",
        show_options = False, 
        link = False
    ):
    
    node = group.nodes.new('GeometryNodeGroup')
    node.node_tree = append(name, link = link)
    
    # if there is an input socket called 'Material', assign it to the base MN material
    # if another material is supplied, use that instead.
    assign_material(node, material=material)
    
    # move and format the node for arranging
    node.location = location
    node.width = width
    node.show_options = show_options
    node.name = name
    node.label = format_node_name(name)
    
    return node

def change_style_node(object, style):
    group = get_mod(object).node_group
    link = group.links.new
    node_style = get_style_node(object)
    socket_from = node_style.inputs[0].links[0].from_socket
    socket_to   = node_style.outputs[0].links[0].to_socket
    new_style = append(styles_mapping[style])
    node_style.node_tree = new_style
    node_style.name = new_style.name
    node_style.label = format_node_name(new_style.name)
    link(socket_from, node_style.inputs[0])
    link(node_style.outputs[0], socket_to)
    assign_material(get_style_node(object))

def create_starting_nodes_starfile(object):
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    node_mod = get_mod(object)
    
    node_name = f"MN_starfile_{object.name}"
    
    # Make sure the aotmic material is loaded
    MN_base_material()
    # create a new GN node group, specific to this particular molecule
    group = new_group(node_name)
    node_mod.node_group = group
    link = group.links.new

     # move the input and output nodes for the group
    node_input = get_input(group)
    node_output = get_output(group)
    node_input.location = [0, 0]
    node_output.location = [700, 0]
    node_star_instances = add_custom(group, 'MN_starfile_instances', [450, 0])
    link(node_star_instances.outputs[0], node_output.inputs[0])
    link(node_input.outputs[0], node_star_instances.inputs[0])
    
    # Need to manually set Image input to 1, otherwise it will be 0 (even though default is 1)
    node_mod['Input_3'] = 1

def create_starting_nodes_density(object, threshold = 0.8, style = 'surface'):
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    mod = get_mod(object)
    node_name = f"MN_density_{object.name}"
    
    # create a new GN node group, specific to this particular molecule
    group = new_group(node_name)
    link = group.links.new
    mod.node_group = group
    
    # move the input and output nodes for the group
    node_input = get_input(group)
    node_input.location = [0, 0]
    node_output = get_output(group)
    node_output.location = [800, 0]
    
    node_density = add_custom(group, styles_mapping[style], [400, 0])
    node_density.inputs['Threshold'].default_value = threshold
    
    link(node_input.outputs[0], node_density.inputs[0])
    link(node_density.outputs[0], node_output.inputs[0])


def create_starting_node_tree(object, coll_frames=None, style="spheres", name=None, set_color=True):
    """
    Create a starting node tree for the inputted object.

    Parameters
    ----------
    object : bpy.types.Object
        Object to create the node tree for.
    coll_frames : bpy.data.collections, optional
        If None, no animation will be created.
        The default is None.
    style : str, optional
        Starting style for the node tree. The default is "spheres".
        Available options are stored as the keys of styles_mapping.
    name : str, optional
        Name of the node tree. If None, a default name will be generated.
        The default is None.
    set_color : bool, optional
        Whether to set up nodes for generating colors in the node tree.
        The default is True.
    """
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    mod = get_mod(object)
    
    if not name:
        name = f"MN_{object.name}"
    
    # create a new GN node group, specific to this particular molecule
    group = new_group(name)
    link = group.links.new
    mod.node_group = group
    
    # move the input and output nodes for the group
    node_input = get_input(group)
    node_output = get_output(group)
    node_input.location = [0, 0]
    node_output.location = [700, 0]
    
    node_style = add_custom(group, styles_mapping[style], [450, 0])
    link(node_style.outputs[0], node_output.inputs[0])
    link(node_input.outputs[0], node_style.inputs[0])
    
    # if requested, setup the nodes for generating colors in the node tree
    if set_color:
        node_color_set = add_custom(group, 'MN_color_set', [200, 0])
        node_color_common = add_custom(group, 'MN_color_common', [-50, -150])
        node_random_color = add_custom(group, 'MN_color_attribute_random', [-300, -150])
        
        link(node_input.outputs['Geometry'], node_color_set.inputs[0])
        link(node_random_color.outputs['Color'], node_color_common.inputs['Carbon'])
        link(node_color_common.outputs[0], node_color_set.inputs['Color'])    
        link(node_color_set.outputs[0], node_style.inputs[0])
        to_animate = node_color_set
    else:
        to_animate = node_input

    # if multiple frames, set up the required nodes for an animation
    if coll_frames:
        node_output.location = [1100, 0]
        node_style.location = [800, 0]
        
        node_animate_frames = add_custom(group, 'MN_animate_frames', [500, 0])
        node_animate = add_custom(group, 'MN_animate_value', [500, -300])
        
        node_animate_frames.inputs['Frames'].default_value = coll_frames
        node_animate.inputs['To Max'].default_value = len(coll_frames.objects) - 1
        
        link(to_animate.outputs[0], node_animate_frames.inputs[0])
        link(node_animate_frames.outputs[0], node_style.inputs[0])
        link(node_animate.outputs[0], node_animate_frames.inputs['Frame'])


def combine_join_geometry(group, node_list, output = 'Geometry', join_offset = 300):
    link = group.links.new
    max_x = max([node.location[0] for node in node_list])
    node_to_instances = group.nodes.new('GeometryNodeJoinGeometry')
    node_to_instances.location = [int(max_x + join_offset), 0]

    for node in reversed(node_list):
        link(node.outputs[output], node_to_instances.inputs['Geometry'])
    return node_to_instances


def split_geometry_to_instances(name, iter_list=('A', 'B', 'C'), attribute='chain_id'):
    """Create a Node to Split Geometry by an Attribute into Instances
    
    Splits the inputted geometry into instances, based on an attribute field. By
    default this field is the `chain_id` but this can be selected for any field.
    Will loop over each item of the list, so a list of arbitrary items that will 
    define how many times to create the required nodes.
    
    """
    group = new_group(name)
    node_input = get_input(group)
    node_output = get_output(group)

    named_att = group.nodes.new('GeometryNodeInputNamedAttribute')
    named_att.location = [-200, -200]
    named_att.data_type = 'INT'
    named_att.inputs[0].default_value = attribute

    link = group.links.new
    list_sep = []

    for i, chain in enumerate(iter_list):
        pos = [i % 10, math.floor(i / 10)]

        node_split = add_custom(group, '.MN_utils_split_instance')
        node_split.location = [int(250 * pos[0]), int(-300 * pos[1])]
        node_split.inputs['Group ID'].default_value = i
        link(named_att.outputs[4], node_split.inputs['Field'])
        link(node_input.outputs['Geometry'], node_split.inputs['Geometry'])
        list_sep.append(node_split)

    node_instance = combine_join_geometry(group, list_sep, 'Instance')
    node_output.location = [int(10 * 250 + 400), 0]
    link(node_instance.outputs[0], node_output.inputs[0])
    return group

def assembly_initialise(mol: bpy.types.Object):
    """
    Setup the required data object and nodes for building an assembly.
    """
    
    transforms = assembly.mesh.array_quaternions_from_dict(mol['biological_assemblies'])
    data_object = assembly.mesh.create_data_object(
        transforms_array=transforms, 
        name = f"data_assembly_{mol.name}"
    )
    tree_assembly = create_assembly_node_tree(
        name = mol.name, 
        iter_list = mol['chain_id_unique'], 
        data_object = data_object
    )
    return tree_assembly

def assembly_insert(mol: bpy.types.Object):
    """
    Given a molecule, setup the required assembly node and insert it into the node tree.
    """
    
    tree_assembly = assembly_initialise(mol)
    group = get_mod(mol).node_group
    node = add_custom(group, tree_assembly.name)
    insert_last_node(get_mod(mol).node_group, node)

def create_assembly_node_tree(name, iter_list, data_object):
    
    node_group_name = f"MN_assembly_{name}"
    group = new_group(name = node_group_name)
    link = group.links.new
    
    n_assemblies = len(np.unique(obj.get_attribute(data_object, 'assembly_id')))
    
    node_group_instances = split_geometry_to_instances(
        name = f".MN_utils_split_{name}", 
        iter_list = iter_list, 
        attribute = 'chain_id'
    )
    
    node_group_assembly_instance = append('.MN_assembly_instance_chains')
    node_instances = add_custom(group, node_group_instances.name, [0, 0])
    node_assembly = add_custom(group, node_group_assembly_instance.name, [200, 0])
    node_assembly.inputs['data_object'].default_value = data_object
    
    out_sockets = outputs(group)
    out_sockets[list(out_sockets)[0]].name = "Instances"
    
    socket_info = (
        {"name" : "Rotation",    "type": "NodeSocketFloat", "min": 0, "max": 1, "default": 1},
        {"name" : "Translation", "type": "NodeSocketFloat", "min": 0, "max": 1, "default": 1},
        {"name" : "assembly_id", "type": "NodeSocketInt", "min": 1, "max": n_assemblies, "default": 1}
    )
    
    for info in socket_info:
        socket = group.interface.items_tree.get(info['name'])
        if not socket:
            socket = group.interface.new_socket(info['name'], in_out='INPUT', socket_type=info['type'])
        socket.default_value = info['default']
        socket.min_value = info['min']
        socket.max_value = info['max']
        
        link(get_input(group).outputs[info['name']], node_assembly.inputs[info['name']])
    
    get_output(group).location = [400, 0]
    link(get_input(group).outputs[0], node_instances.inputs[0])
    link(node_instances.outputs[0], node_assembly.inputs[0])
    link(node_assembly.outputs[0], get_output(group).inputs[0])
    
    return group


def add_inverse_selection(group):
    output = get_output(group)
    if 'Inverted' not in output.inputs.keys():
        group.interface.new_socket('Inverted', in_out='OUTPUT', socket_type='NodeSocketBool')
    
    loc = output.location
    bool_math = group.nodes.new("FunctionNodeBooleanMath")
    bool_math.location = [loc[0], -100]
    bool_math.operation = "NOT"
    
    group.links.new(output.inputs['Selection'].links[0].from_socket, bool_math.inputs[0])
    group.links.new(bool_math.outputs[0], output.inputs['Inverted'])

def chain_selection(name, input_list, attribute = 'chain_id', starting_value = 0, label_prefix = ""):
    """
    Given a an input_list, will create a node which takes an Integer input, 
    and has a boolean tick box for each item in the input list. The outputs will
    be the resulting selection and the inversion of the selection.
    Can contain a prefix for the resulting labels. Mostly used for constructing 
    chain selections when required for specific proteins.
    """
    # just reutn the group name early if it already exists
    group = bpy.data.node_groups.get(name)
    if group:
        return group
    
    group = new_group(name, geometry=False)
    link = group.links.new
    # create a named attribute node that gets the chain_number attribute
    # and use this for the selection algebra that happens later on
    node_attribute = group.nodes.new("GeometryNodeInputNamedAttribute")
    node_attribute.data_type = 'INT'
    node_attribute.location = [-200, 200]
    node_attribute.inputs[0].default_value = attribute
    node_attribute.outputs.get('Attribute')
    
    # distance horizontally to space all of the created nodes
    node_sep_dis = 180
    previous_node = None
    att_output = get_output_type(node_attribute, 'INT')
    
    for i, chain_name in enumerate(input_list):
        group.interface.new_socket(str(label_prefix) + str(chain_name), in_out='INPUT', socket_type='NodeSocketBool')
        current_node = group.nodes.new("GeometryNodeGroup")
        current_node.node_tree = append('.MN_utils_bool_chain')
        current_node.location = [i * node_sep_dis, 200]
        current_node.inputs["number_matched"].default_value = i + starting_value

        # link from the the named attribute node chain_number into the other inputs
        link(get_input(group).outputs[i], current_node.inputs["bool_include"])
        if previous_node:
            link(previous_node.outputs['number_chain_out'], current_node.inputs['number_chain_in'])
            link(previous_node.outputs['bool_chain_out'], current_node.inputs['bool_chain_in'])
        else:
            link(att_output, current_node.inputs['number_chain_in'])
        previous_node = current_node

    group.interface.new_socket('Selection', in_out='OUTPUT', socket_type='NodeSocketBool')
    group_out = get_output(group)
    group_out.location = [len(input_list) * node_sep_dis, 200]
    link(current_node.outputs['bool_chain_out'], group_out.inputs['Selection'])
    add_inverse_selection(group)

    return group

def chain_color(name, input_list, label_prefix = "Chain ", field = "chain_id", starting_value = 0):
    """
    Given the input list of chain names, will create a node group which uses
    the chain_id named attribute to manually set the colours for each of the chains.
    """
    
    group = new_group(name, geometry=False)
    link = group.links.new
    
    # create a named attribute node that gets the chain_number attribute
    # and use this for the selection algebra that happens later on
    node_att = group.nodes.new("GeometryNodeInputNamedAttribute")
    node_att.data_type = 'INT'
    node_att.location = [-200, 400]
    node_att.inputs[0].default_value = field
    node_att.outputs.get('Attribute')
    
    node_input = get_input(group)
    node_output = get_output(group)
    
    # shortcut for creating new nodes
    new_node = group.nodes.new
    # distance horizontally to space all of the created nodes
    node_sep_dis = 180
    att_output = get_output_type(node_att, 'INT')
    node_color_previous = None
    
    for i, chain_name in enumerate(input_list):
        i += starting_value
        offset = i * node_sep_dis
        current_chain = f"{label_prefix}{chain_name}"
        
        # node compare inputs 2 & 3
        node_compare = new_node('FunctionNodeCompare')
        node_compare.data_type = 'INT'
        node_compare.location = [offset, 100]
        node_compare.operation = 'EQUAL'
        
        node_compare.inputs[3].default_value = i
        
        # link the named attribute to the compare
        link(att_output, node_compare.inputs[2])
        
        node_color = new_node('GeometryNodeSwitch')
        node_color.input_type = 'RGBA'
        node_color.location = [offset, -100]
        
        # create an input for this chain
        group.interface.new_socket(current_chain, in_out='INPUT', socket_type='NodeSocketColor').default_value = color.random_rgb(i)
        # switch input colours 10 and 11
        link(node_input.outputs[current_chain], node_color.inputs[11])
        link(node_compare.outputs['Result'], node_color.inputs['Switch'])
        
        
        if node_color_previous:
            link(node_color_previous.outputs[4], node_color.inputs[10])
        node_color_previous = node_color
    
    group.interface.new_socket('Color', in_out='OUTPUT', socket_type='NodeSocketColor')
    link(node_color.outputs[4], node_output.inputs['Color'])
    
    return group

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
    
    # create the custom node group data block, where everything will go
    # also create the required group node input and position it
    residue_id_group = bpy.data.node_groups.new(node_name, "GeometryNodeTree")
    node_input = residue_id_group.nodes.new("NodeGroupInput")
    node_input.location = [0, node_sep_dis * len(sub_list)/2]
    
    group_link = residue_id_group.links.new
    new_node = residue_id_group.nodes.new
    
    for residue_id_index,residue_id in enumerate(sub_list):
        
        # add an new node of Select Res ID or MN_sek_res_id_range
        current_node = new_node("GeometryNodeGroup")

        # add an bool_math block 
        bool_math = new_node("FunctionNodeBooleanMath")
        bool_math.location = [400,(residue_id_index+1) * node_sep_dis]
        bool_math.operation = "OR"

        if '-' in residue_id:
            # set two new inputs
            current_node.node_tree = append('MN_select_res_id_range')
            [resid_start, resid_end] = residue_id.split('-')[:2]
            socket_1 = residue_id_group.interface.new_socket('res_id: Min', in_out='INPUT', socket_type='NodeSocketInt')
            socket_1.default_value = int(resid_start)
            socket_2 = residue_id_group.interface.new_socket('res_id: Max', in_out='INPUT', socket_type='NodeSocketInt')
            socket_2.default_value = int(resid_end)
            
            # a residue range
            group_link(node_input.outputs[socket_1.identifier], current_node.inputs[0])
            group_link(node_input.outputs[socket_2.identifier], current_node.inputs[1])
        else:
            # create a node
            current_node.node_tree = append('MN_select_res_id_single')
            socket = residue_id_group.interface.new_socket('res_id', in_out='INPUT', socket_type='NodeSocketInt')
            socket.default_value  = int(residue_id)
            group_link(node_input.outputs[socket.identifier], current_node.inputs[0])
        
        # set the coordinates
        current_node.location = [200,(residue_id_index+1) * node_sep_dis]
        prev = None
        if not prev:
            # link the first residue selection to the first input of its OR block
            group_link(current_node.outputs['Selection'],bool_math.inputs[0])
        else:
            # if it is not the first residue selection, link the output to the previous or block
            group_link(current_node.outputs['Selection'], prev.inputs[1])
        
            # link the ouput of previous OR block to the current OR block
            group_link(prev.outputs[0], bool_math.inputs[0])
        prev = bool_math

    # add a output block
    residue_id_group_out = new_node("NodeGroupOutput")
    residue_id_group_out.location = [800,(residue_id_index + 1) / 2 * node_sep_dis]
    residue_id_group.interface.new_socket('Selection', in_out='OUTPUT', socket_type='NodeSocketBool')
    residue_id_group.interface.new_socket('Inverted', in_out='OUTPUT', socket_type='NodeSocketBool')
    group_link(prev.outputs[0], residue_id_group_out.inputs['Selection'])
    invert_bool_math = new_node("FunctionNodeBooleanMath")
    invert_bool_math.location = [600,(residue_id_index+1)/ 3 * 2 * node_sep_dis]
    invert_bool_math.operation = "NOT"
    group_link(prev.outputs[0], invert_bool_math.inputs[0])
    group_link(invert_bool_math.outputs[0], residue_id_group_out.inputs['Inverted'])
    return residue_id_group