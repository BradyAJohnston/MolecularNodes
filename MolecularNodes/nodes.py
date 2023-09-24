import bpy
import os
from . import pkg
import math
from . import obj
import numpy as np
from . import color
import warnings

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
    'atoms': 'MN_style_atoms',
    'vdw': 'MN_style_atoms',
    'sphere': 'MN_style_atoms',
    'cartoon': 'MN_style_cartoon',
    'ribbon': 'MN_style_ribbon_protein',
    'ball_and_stick': 'MN_style_ball_and_stick',
    'ball+stick': 'MN_style_ball_and_stick',
}


mn_data_file = os.path.join(pkg.ADDON_DIR, 'assets', 'MN_data_file.blend')

def append(node_name, link = False):
    node = bpy.data.node_groups.get(node_name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if not node or link:
            bpy.ops.wm.append(
                directory = os.path.join(mn_data_file, 'NodeTree'), 
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
            directory = os.path.join(mn_data_file, 'Material'), 
            filename = 'MN_atomic_material', 
            link = False
        )
    
    return bpy.data.materials[mat_name]

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

def add_node(node_name, label: str = '', show_options = False):
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
    if label != '':
        node.label = label
    
    # if added node has a 'Material' input, set it to the default MN material
    input_mat = bpy.context.active_node.inputs.get('Material')
    if input_mat:
        input_mat.default_value = MN_base_material()

def add_custom_node_group(parent_group, 
                          node_name, 
                          location = [0,0], 
                          width = 200, 
                          show_options = False, 
                          link = False
                          ):
    
    append(node_name, link=link)
    
    node = parent_group.node_group.nodes.new('GeometryNodeGroup')
    node.node_tree = bpy.data.node_groups[node_name]
    
    node.location = location
    node.width = 200 # unsure if width will work TODO check
    node.show_options = show_options
    node.name = node_name
    
    return node

def add_custom_node_group_to_node(parent_group, 
                                  node_name, 
                                  location = [0,0], 
                                  width = 200, 
                                  show_options = False, 
                                  link = False
                                  ):
    
    append(node_name, link = link)
    
    node = parent_group.nodes.new('GeometryNodeGroup')
    node.node_tree = bpy.data.node_groups[node_name]
    
    node.location = location
    node.width = 200 # unsure if width will work TODO check
    node.show_options = show_options
    node.name = node_name
    
    return node

def create_starting_nodes_starfile(obj):
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    node_mod = obj.modifiers.get('MolecularNodes')
    if not node_mod:
        node_mod = obj.modifiers.new("MolecularNodes", "NODES")
    obj.modifiers.active = node_mod
    
    node_name = f"MN_starfile_{obj.name}"
    
    # if node tree already exists by this name, set it and return it
    node_group = bpy.data.node_groups.get(node_name)
    if node_group:
        node_mod.node_group = node_group
        return node_group
    
    
    # create a new GN node group, specific to this particular molecule
    node_group = gn_new_group_empty(node_name)
    
    # create a new GN node group, specific to this particular molecule
    node_group = gn_new_group_empty(node_name)
    node_mod.node_group = node_group
    node_group.inputs.new("NodeSocketObject", "Molecule")
    node_group.inputs.new("NodeSocketInt", "Image")
    node_group.inputs["Image"].default_value = 1
    node_group.inputs["Image"].min_value = 1
    node_group.inputs.new("NodeSocketBool", "Simplify")
    # move the input and output nodes for the group
    node_input = node_mod.node_group.nodes[bpy.app.translations.pgettext_data("Group Input",)]
    node_input.location = [0, 0]
    node_output = node_mod.node_group.nodes[bpy.app.translations.pgettext_data("Group Output",)]
    node_output.location = [900, 0]

    node_delete = node_group.nodes.new("GeometryNodeDeleteGeometry")
    node_delete.location = [500, 0]

    node_geom_to_instance = node_group.nodes.new("GeometryNodeInstanceOnPoints")
    node_geom_to_instance.location = [675, 0]

    node_get_imageid = node_group.nodes.new("GeometryNodeInputNamedAttribute")
    node_get_imageid.location = [0, 200]
    node_get_imageid.inputs['Name'].default_value = "MOLImageId"
    node_get_imageid.data_type = "INT"

    node_subtract = node_group.nodes.new("ShaderNodeMath")
    node_subtract.location = [160, 200]
    node_subtract.operation = "SUBTRACT"
    node_subtract.inputs[1].default_value = 1
    node_subtract.inputs[0].default_value = 1


    node_compare = node_group.nodes.new("FunctionNodeCompare")
    node_compare.location = [320, 200]
    node_compare.operation = "NOT_EQUAL"
    node_compare.data_type = "INT"

    node_object_info = node_group.nodes.new("GeometryNodeObjectInfo")
    node_object_info.location = [200, -200]

    node_get_rotation = node_group.nodes.new("GeometryNodeInputNamedAttribute")
    node_get_rotation.location = [450, -200]
    node_get_rotation.inputs['Name'].default_value = "MOLRotation"
    node_get_rotation.data_type = "FLOAT_VECTOR"

    node_get_id = node_group.nodes.new("GeometryNodeInputID")
    node_get_id.location = [0, -200]

    node_statistics = node_group.nodes.new("GeometryNodeAttributeStatistic")
    node_statistics.location = [200, -400]

    node_compare_maxid = node_group.nodes.new("FunctionNodeCompare")
    node_compare_maxid.location = [400, -400]
    node_compare_maxid.operation = "EQUAL"

    node_bool_math = node_group.nodes.new("FunctionNodeBooleanMath")
    node_bool_math.location = [600, -400]
    node_bool_math.operation = "OR"

    node_switch = node_group.nodes.new("GeometryNodeSwitch")
    node_switch.location = [800, -400]

    node_cone = node_group.nodes.new("GeometryNodeMeshCone")
    node_cone.location = [1000, -400]

    link = node_group.links.new

    link(node_input.outputs[0], node_delete.inputs[0])
    link(node_delete.outputs[0], node_geom_to_instance.inputs[0])
    link(node_geom_to_instance.outputs[0], node_output.inputs[0])

    link(node_input.outputs[1], node_object_info.inputs[0])
    link(node_input.outputs[2], node_subtract.inputs[0])
    link(node_input.outputs[3], node_bool_math.inputs[0])

    link(node_subtract.outputs[0], node_compare.inputs[2])
    link(node_get_imageid.outputs[4], node_compare.inputs[3])
    link(node_compare.outputs[0], node_delete.inputs[1])
    link(node_statistics.outputs[4], node_compare_maxid.inputs[0])
    link(node_compare_maxid.outputs[0], node_bool_math.inputs[1])
    link(node_get_id.outputs[0], node_statistics.inputs[2])
    link(node_object_info.outputs["Geometry"], node_statistics.inputs[0])
    link(node_bool_math.outputs[0], node_switch.inputs[1])
    link(node_object_info.outputs["Geometry"], node_switch.inputs[14])
    link(node_cone.outputs[0], node_switch.inputs[15])
    link(node_switch.outputs[6],     node_geom_to_instance.inputs["Instance"])
    link(node_get_rotation.outputs[0], node_geom_to_instance.inputs["Rotation"])


    # Need to manually set Image input to 1, otherwise it will be 0 (even though default is 1)
    node_mod['Input_3'] = 1

def create_starting_nodes_density(obj, threshold = 0.8):
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    node_mod = obj.modifiers.get('MolecularNodes')
    if not node_mod:
        node_mod = obj.modifiers.new("MolecularNodes", "NODES")
    obj.modifiers.active = node_mod
    node_name = f"MN_density_{obj.name}"
    
    # if node tree already exists by this name, set it and return it
    node_group = bpy.data.node_groups.get(node_name)
    if node_group:
        node_mod.node_group = node_group
        return node_group
    
    
    # create a new GN node group, specific to this particular molecule
    node_group = gn_new_group_empty(node_name)
    node_mod.node_group = node_group
    # move the input and output nodes for the group
    node_input = node_mod.node_group.nodes[bpy.app.translations.pgettext_data("Group Input",)]
    node_input.location = [0, 0]
    node_output = node_mod.node_group.nodes[bpy.app.translations.pgettext_data("Group Output",)]
    node_output.location = [800, 0]
    
    node_density = add_custom_node_group(node_mod, 'MN_density_style_surface', [400, 0])
    node_density.inputs['Material'].default_value = MN_base_material()
    node_density.inputs['Density Threshold'].default_value = threshold
    
    
    link = node_group.links.new
    link(
        node_input.outputs[0], 
        node_density.inputs[0]
    )
    link(
        node_density.outputs[0], 
        node_output.inputs[0]
    )
    
    

def create_starting_node_tree(obj,
                              coll_frames=None,
                              starting_style="atoms"):
    """
    Create a starting node tree for the inputted object.

    Parameters
    ----------
    obj : bpy.types.Object
        Object to create the node tree for.
    coll_frames : bpy.data.collections, optional
        If None, no animation will be created.
        The default is None.
    starting_style : str, optional
        Starting style for the node tree. The default is "atoms".
        Available options are stored as the keys of styles_mapping
    """
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    node_mod = obj.modifiers.get('MolecularNodes')
    if not node_mod:
        node_mod = obj.modifiers.new("MolecularNodes", "NODES")
    obj.modifiers.active = node_mod
    
    
    name = f"MN_{obj.name}"
    # if node group of this name already exists, set that node group
    # and return it without making any changes
    node_group = bpy.data.node_groups.get(name)
    if node_group:
        node_mod.node_group = node_group
        return node_group
    
    # create a new GN node group, specific to this particular molecule
    node_group = gn_new_group_empty(name)
    node_mod.node_group = node_group
    
    # move the input and output nodes for the group
    node_input = node_mod.node_group.nodes[bpy.app.translations.pgettext_data("Group Input",)]
    node_input.location = [0, 0]
    node_output = node_mod.node_group.nodes[bpy.app.translations.pgettext_data("Group Output",)]
    node_output.location = [700, 0]
    
    # node_properties = add_custom_node_group(node_group, 'MN_prop_setup', [0, 0])
    node_color_set = add_custom_node_group(node_mod, 'MN_color_set', [200, 0])
    node_color_common = add_custom_node_group(node_mod, 'MN_color_common', [-50, -150])
    
    
    
    node_random_color = add_custom_node_group(node_mod, 'MN_color_attribute_random', [-300, -150])
    
    
    # create the links between the the nodes that have been established
    link = node_group.links.new
    link(node_input.outputs['Geometry'], node_color_set.inputs[0])
    link(node_color_set.outputs[0], node_output.inputs['Geometry'])
    link(node_random_color.outputs['Color'], node_color_common.inputs['Carbon'])
    link(node_color_common.outputs[0], node_color_set.inputs['Color'])

    
    node_style = add_custom_node_group(node_mod,
                                       styles_mapping[starting_style],
                                       location = [450, 0])
    link(node_color_set.outputs['Atoms'], node_style.inputs['Atoms'])
    link(node_style.outputs[0], node_output.inputs['Geometry'])
    node_style.inputs['Material'].default_value = MN_base_material()

    # if multiple frames, set up the required nodes for an animation
    if coll_frames:
        node_output.location = [1100, 0]
        node_style.location = [800, 0]
        
        node_animate_frames = add_custom_node_group_to_node(node_group, 'MN_animate_frames', [500, 0])
        node_animate_frames.inputs['Frames'].default_value = coll_frames
        
        # node_animate_frames.inputs['Absolute Frame Position'].default_value = True
        
        node_animate = add_custom_node_group_to_node(node_group, 'MN_animate_value', [500, -300])
        link(node_color_set.outputs['Atoms'], node_animate_frames.inputs['Atoms'])
        link(node_animate_frames.outputs['Atoms'], node_style.inputs['Atoms'])
        link(node_animate.outputs[0], node_animate_frames.inputs['Animate 0..1'])


def split_geometry_to_instances(name, iter_list=('A', 'B', 'C'), attribute='chain_id'):
    """Create a Node to Split Geometry by an Attribute into Instances
    
    Splits the inputted geometry into instances, based on an attribute field. By
    default this field is the `chain_id` but this can be selected for any field.
    Will loop over each item of the list, so a list of arbitrary items that will 
    define how many times to create the required nodes.
    
    """
    node_group = bpy.data.node_groups.get(name)
    if node_group:
        return node_group

    node_group = gn_new_group_empty(name)

    node_input = node_group.nodes[bpy.app.translations.pgettext_data('Group Input')]
    node_output = node_group.nodes[bpy.app.translations.pgettext_data('Group Output')]

    add_node = node_group.nodes.new
    named_att = add_node('GeometryNodeInputNamedAttribute')
    named_att.location = [-200, -200]
    named_att.data_type = 'INT'
    named_att.inputs[0].default_value = attribute

    link = node_group.links.new
    list_sep = []

    for i, chain in enumerate(iter_list):

        pos = [i % 10, math.floor(i / 10)]

        node_split = add_custom_node_group_to_node(node_group, '.MN_utils_split_instance')
        node_split.location = [int(250 * pos[0]), int(-300 * pos[1])]
        node_split.inputs['Group ID'].default_value = i


        link(named_att.outputs[4], node_split.inputs['Field'])
        link(node_input.outputs['Geometry'], node_split.inputs['Geometry'])

        list_sep.append(node_split)

    node_instance = nodes_to_geometry(node_group, list_sep, 'Instance')

    node_output.location = [int(10 * 250 + 400), 0]
    link(node_instance.outputs[0], node_output.inputs[0])

    return node_group

def nodes_to_geometry(this_group, node_list, output = 'Geometry', join_offset = 300):
    link = this_group.links.new
    max_x = max([node.location[0] for node in node_list])
    node_to_instances = this_group.nodes.new('GeometryNodeJoinGeometry')
    node_to_instances.location = [int(max_x + join_offset), 0]

    for node in reversed(node_list):
        link(
            node.outputs[output], 
            node_to_instances.inputs['Geometry']
        )

    return node_to_instances

def create_assembly_node_tree(name, iter_list, data_object):
    
    node_group_name = f"MN_assembly_{name}"
    group = bpy.data.node_groups.get(node_group_name)
    if group:
        return group
    
    group = gn_new_group_empty(name = node_group_name)
    
    n_assemblies = len(np.unique(obj.get_attribute(data_object, 'assembly_id')))
    
    node_group_instances = split_geometry_to_instances(
        name = f"MN_utils_split_{name}", 
        iter_list = iter_list, 
        attribute = 'chain_id'
    )
    
    node_group_assembly_instance = append('.MN_assembly_instance_chains')
    
    def new_node_group(name, location = [0, 0]):
        node = group.nodes.new("GeometryNodeGroup")
        node.node_tree = bpy.data.node_groups[name]
        node.location = location
        return node
    
    link = group.links.new
    
    node_instances = new_node_group(node_group_instances.name, [0, 0])
    node_assembly = new_node_group(node_group_assembly_instance.name, [200, 0])
    node_assembly.inputs['data_object'].default_value = data_object
    
    group.outputs[0].name = "Assembly Instances"
    
    socket_info = (
        {"name" : "Rotation",    "type": "NodeSocketFloat", "min": 0, "max": 1, "default": 1},
        {"name" : "Translation", "type": "NodeSocketFloat", "min": 0, "max": 1, "default": 1},
        {"name" : "assembly_id", "type": "NodeSocketInt", "min": 1, "max": n_assemblies, "default": 1}
    )
    
    for socket in socket_info:
        new_socket = group.inputs.get(socket['name'])
        if not new_socket:
            new_socket = group.inputs.new(socket['type'], socket['name'])
        new_socket.default_value = socket['default']
        new_socket.min_value = socket['min']
        new_socket.max_value = socket['max']
        
        link(
            group.nodes['Group Input'].outputs[socket['name']], 
            node_assembly.inputs[socket['name']]
        )
    
    group.nodes['Group Output'].location = [400, 0]
    
    link(
        group.nodes['Group Input'].outputs[0], 
        node_instances.inputs[0]
    )
    link(
        node_instances.outputs[0], 
        node_assembly.inputs[0]
    )
    link(
        node_assembly.outputs[0], 
        group.nodes['Group Output'].inputs[0]
    )
    
    return group

def create_custom_surface(name, n_chains):
    # if a group of this name already exists, just return that instead of making a new one
    group = bpy.data.node_groups.get(name)
    if group:
        return group
    
    # get the node to create a loop from
    looping_node = append('MN_style_surface')
    
    
    # create new empty data block
    group = bpy.data.node_groups.new(name, 'GeometryNodeTree')
    
    # loop over the inputs and create an input for each
    for i in looping_node.inputs.values():
        group_input = group.inputs.new(socket_types.get(i.type), i.name)
        try:
            group_input.default_value = i.default_value
        except AttributeError:
            pass
    
    group.inputs['Selection'].hide_value = True
    
    group.outputs.new(socket_types.get('GEOMETRY'), 'Surface Geometry')
    group.outputs.new(socket_types.get('GEOMETRY'), 'Surface Instances')
    
    # add in the inputs and theo outputs inside of the node
    node_input = group.nodes.new('NodeGroupInput')
    node_input.location = [-300, 0]
    node_output = group.nodes.new('NodeGroupOutput')
    node_output.location = [800, 0]
    
    link = group.links.new
    
    node_input = group.nodes[bpy.app.translations.pgettext_data("Group Input",)]
    # node_output = group.nodes[bpy.app.translations.pgettext_data("Group Output",)]
    
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
        
        for i in node_surface_single.inputs.values():
            if i.type != 'GEOMETRY':
                link(node_input.outputs[i.name], i)
        
        
        list_node_surface.append(node_surface_single)
    
    # create join geometry, and link the nodes in reverse order
    node_join_geometry = group.nodes.new('GeometryNodeJoinGeometry')
    node_join_geometry.location = [500, 0]
    
    node_geom_to_instance = group.nodes.new('GeometryNodeGeometryToInstance')
    node_geom_to_instance.location = [500, -600]
    
    for n in reversed(list_node_surface): 
        link(n.outputs[0], node_join_geometry.inputs['Geometry'])
        link(n.outputs[0], node_geom_to_instance.inputs['Geometry'])
    
    # link the joined nodes to the outputs
    link(node_join_geometry.outputs['Geometry'], node_output.inputs[0])
    link(node_geom_to_instance.outputs['Instances'], node_output.inputs['Surface Instances'])
    
    return group

def chain_selection(node_name, input_list, attribute = 'chain_id', starting_value = 0, label_prefix = ""):
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
        current_node.node_tree = append('.MN_utils_bool_chain')
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

def chain_color(node_name, input_list, label_prefix = "Chain ", field = "chain_id", starting_value = 0):
    """
    Given the input list of chain names, will create a node group which uses
    the chain_id named attribute to manually set the colours for each of the chains.
    """
    
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
    chain_number_node.inputs[0].default_value = field
    chain_number_node.outputs.get('Attribute')
    
    # shortcut for creating new nodes
    new_node = chain_group.nodes.new
    # distance horizontally to space all of the created nodes
    node_sep_dis = 180
    counter = starting_value
    
    for chain_name in input_list:
        offset = counter * node_sep_dis
        current_chain = f"{label_prefix}{chain_name}"
        
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
        chain_group.inputs[current_chain].default_value = color.random_rgb()
        # switch input colours 10 and 11
        link(node_input.outputs[current_chain], node_color.inputs[11])
        link(node_compare.outputs['Result'], node_color.inputs['Switch'])
        
        
        if counter > starting_value:
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
            residue_id_group.inputs.new("NodeSocketInt",'res_id: Min').default_value = int(resid_start)
            residue_id_group.inputs.new("NodeSocketInt",'res_id: Max').default_value = int(resid_end)
        else:
            # set a new input and set the resid
            residue_id_group.inputs.new("NodeSocketInt",'res_id').default_value = int(residue_id)
        
    # set a counter for Select Res ID* nodes
    counter=0
    for residue_id_index,residue_id in enumerate(sub_list):
        
        # add an new node of Select Res ID or MN_sek_res_id_range
        current_node = new_node("GeometryNodeGroup")

        # add an bool_math block 
        bool_math = new_node("FunctionNodeBooleanMath")
        bool_math.location = [400,(residue_id_index+1) * node_sep_dis]
        bool_math.operation = "OR"

        if '-' in residue_id:
            # a residue range
            current_node.node_tree = append('MN_select_res_id_range')
            
            group_link(residue_id_group_in.outputs[counter], current_node.inputs[0])
            counter+=1
            
            group_link(residue_id_group_in.outputs[counter], current_node.inputs[1])
            
        else:
            # create a node
            current_node.node_tree = append('MN_select_res_id_single')
            # link the input of MN_select_res_ID*
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