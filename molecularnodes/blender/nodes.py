import bpy
import os
import numpy as np
import math
import warnings
import itertools
from .. import utils
from .. import color
from .. import pkg
from ..blender import obj

socket_types = {
    'BOOLEAN': 'NodeSocketBool',
    'GEOMETRY': 'NodeSocketGeometry',
    'INT': 'NodeSocketInt',
    'MATERIAL': 'NodeSocketMaterial',
    'VECTOR': 'NodeSocketVector',
    'STRING': 'NodeSocketString',
    'VALUE': 'NodeSocketFloat',
    'COLLECTION': 'NodeSocketCollection',
    'TEXTURE': 'NodeSocketTexture',
    'COLOR': 'NodeSocketColor',
    'RGBA': 'NodeSocketColor',
    'IMAGE': 'NodeSocketImage'
}

# current implemented representations
styles_mapping = {
    "presets": "MN_style_presets",
    'preset_1': "MN_style_presets",
    'preset_2': "MN_style_presets",
    'preset_3': "MN_style_presets",
    'preset_4': "MN_style_presets",
    'atoms': 'MN_style_spheres',
    'spheres': 'MN_style_spheres',
    'vdw': 'MN_style_spheres',
    'sphere': 'MN_style_spheres',
    'cartoon': 'MN_style_cartoon',
    'ribbon': 'MN_style_ribbon',
    'surface': 'MN_style_surface',
    'ball_and_stick': 'MN_style_ball_and_stick',
    'ball+stick': 'MN_style_ball_and_stick',
    'oxdna': 'MN_oxdna_style_ribbon',
    "density_surface": "MN_density_style_surface",
    "density_wire": "MN_density_style_wire"
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
    name="Style",
    description="Default style for importing molecules.",
    items=STYLE_ITEMS,
    default='spheres'
)


MN_DATA_FILE = os.path.join(pkg.ADDON_DIR, 'assets', 'MN_data_file.blend')


class NodeGroupCreationError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


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


def get_output_type(node, type="INT"):
    for output in node.outputs:
        if output.type == type:
            return output


def set_selection(group, node, selection):
    pos = node.location
    pos = [pos[0] - 200, pos[1] - 200]
    selection.location = pos
    group.links.new(selection.outputs[0], node.inputs['Selection'])

    return selection


def create_debug_group(name='MolecularNodesDebugGroup'):
    group = new_group(name=name, fallback=False)
    info = group.nodes.new('GeometryNodeObjectInfo')
    group.links.new(info.outputs['Geometry'],
                    group.nodes['Group Output'].inputs[0])
    return group


def add_selection(group, sel_name, input_list, field='chain_id'):
    style = style_node(group)
    sel_node = add_custom(group, custom_iswitch(
        name='selection',
        iter_list=input_list,
        field=field,
        dtype='BOOLEAN'
    ).name)

    set_selection(group, style, sel_node)
    return sel_node


def get_output(group):
    return group.nodes[bpy.app.translations.pgettext_data("Group Output",)]


def get_input(group):
    return group.nodes[bpy.app.translations.pgettext_data("Group Input",)]


def get_mod(object, name='MolecularNodes'):
    node_mod = object.modifiers.get(name)
    if not node_mod:
        node_mod = object.modifiers.new(name, "NODES")
    object.modifiers.active = node_mod
    return node_mod


def format_node_name(name):
    "Formats a node's name for nicer printing."
    return name.strip("MN_").replace("_", " ").title().replace('Dna', 'DNA').replace('Topo ', 'Topology ')


def get_nodes_last_output(group):
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


def annotation_instances_node(group):
    prev = previous_node(get_output(group))
    is_annotation_instances_node = ("MN_annotation_instances" in prev.name)
    while not is_annotation_instances_node:
        prev = previous_node(prev)
        is_annotation_instances_node = ("MN_annotation_instances" in prev.name)
    return prev


def get_annotation_instances_node(object):
    "Walk back through the primary node connections until you find the first style node"
    group = object.modifiers['MolecularNodes'].node_group
    return annotation_instances_node(group)


def get_color_node(object):
    "Walk back through the primary node connections until you find the first style node"
    group = object.modifiers['MolecularNodes'].node_group
    for node in group.nodes:
        if node.name == "MN_color_attribute_random":
            return node


def insert_last_node(group, node, link_input=True):
    last, output = get_nodes_last_output(group)
    link = group.links.new
    location = output.location
    output.location = [location[0] + 300, location[1]]
    node.location = [location[0] - 300, location[1]]
    if link_input:
        link(last.outputs[0], node.inputs[0])
    link(node.outputs[0], output.inputs[0])


def realize_instances(obj):
    group = obj.modifiers['MolecularNodes'].node_group
    realize = group.nodes.new('GeometryNodeRealizeInstances')
    insert_last_node(group, realize)


def append(node_name, link=False):
    node = bpy.data.node_groups.get(node_name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if not node or link:
            bpy.ops.wm.append(
                'EXEC_DEFAULT',
                directory=os.path.join(MN_DATA_FILE, 'NodeTree'),
                filename=node_name,
                link=link,
                use_recursive=True
            )
    node = bpy.data.node_groups.get(node_name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if not node or link:
            node_name_components = node_name.split('_')
            if node_name_components[0] == 'MN':
                data_file = MN_DATA_FILE[:-6] + '_' + \
                    node_name_components[1] + '.blend'
                bpy.ops.wm.append(
                    'EXEC_DEFAULT',
                    directory=os.path.join(data_file, 'NodeTree'),
                    filename=node_name,
                    link=link,
                    use_recursive=True
                )
    return bpy.data.node_groups[node_name]


def material_default():
    """
    Append MN Default to the .blend file it it doesn't already exist,
    and return that material.
    """

    mat_name = 'MN Default'
    mat = bpy.data.materials.get(mat_name)

    if not mat:
        print('appending material')
        bpy.ops.wm.append(
            directory=os.path.join(MN_DATA_FILE, 'Material'),
            filename='MN Default',
            link=False
        )

    return bpy.data.materials[mat_name]


def MN_micrograph_material():
    """
    Append MN_micrograph_material to the .blend file it it doesn't already exist,
    and return that material.
    """

    mat_name = 'MN_micrograph_material'

    return bpy.data.materials[mat_name]


def new_group(name="Geometry Nodes", geometry=True, fallback=True):
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
        group.interface.new_socket(
            'Geometry', in_out='INPUT', socket_type='NodeSocketGeometry')
        group.interface.new_socket(
            'Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
        group.links.new(output_node.inputs[0], input_node.outputs[0])
    return group


def assign_material(node, material='default'):
    material_socket = node.inputs.get('Material')
    if material_socket:
        if not material:
            pass
        elif material == "default":
            material_socket.default_value = material_default()
        else:
            material_socket.default_value = material


def add_node(node_name, label: str = '', show_options=False, material="default"):
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
    location=[0, 0],
    width=200,
    material="default",
    show_options=False,
    link=False
):

    node = group.nodes.new('GeometryNodeGroup')
    node.node_tree = append(name, link=link)

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
    # get the node group that we are working on, to change the specific style node
    group = get_mod(object).node_group
    link = group.links.new
    node_style = get_style_node(object)

    # capture input and output links, so we can rebuild the links based on name
    # and the sockets they were connected to
    # as we collect them, remove the links so they aren't automatically connected
    # when we change the node_tree for the group
    input_links = []
    output_links = []
    for input in node_style.inputs:
        for input_link in input.links:
            input_links.append((input_link.from_socket, input.name))
            group.links.remove(input_link)

    for output in node_style.outputs:
        for output_link in output.links:
            output_links.append((output.name, output_link.to_socket))
            group.links.remove(output_link)

    try:
        material = node_style.inputs['Material'].default_value
    except KeyError:
        material = None
    # append the new node tree, and swap out the tree that is used for the group
    group_new = append(styles_mapping[style])
    node_style.node_tree = group_new
    # do some formatting and cleanup the node name for easier reading
    node_style.name = group_new.name
    node_style.label = format_node_name(node_style.name)

    # rebuild the links based on names of the sockets, not their identifiers
    for input_link in input_links:
        link(input_link[0], node_style.inputs[input_link[1]])
    for output_link in output_links:
        link(node_style.outputs[output_link[0]], output_link[1])

    if material:
        try:
            node_style.inputs['Material'].default_value = material
        except KeyError:
            # the new node doesn't contain a material slot
            pass


def create_starting_nodes_annotation_instances(object, n_images=1):
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    node_mod = get_mod(object)

    node_name = f"MN_annotation_instances_{object.name}"

    # Make sure the aotmic material is loaded
    material_default()
    # create a new GN node group, specific to this particular molecule
    group = new_group(node_name)
    node_mod.node_group = group
    link = group.links.new

    # move the input and output nodes for the group
    node_input = get_input(group)
    node_output = get_output(group)
    node_input.location = [0, 0]
    node_output.location = [700, 0]
    node_annotation_instances = add_custom(group, 'MN_annotation_instances', [450, 0])
    link(node_annotation_instances.outputs[0], node_output.inputs[0])
    link(node_input.outputs[0], node_annotation_instances.inputs[0])

    # Need to manually set Image input to 1, otherwise it will be 0 (even though default is 1)
    node_mod['Input_3'] = 1


def create_starting_nodes_density(object, threshold=0.8, style='density_surface'):
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    mod = get_mod(object)
    node_name = f"MN_density_{object.name}"

    # create a new GN node group, specific to this particular molecule
    group = new_group(node_name, fallback=False)
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
    mod = get_mod(object)
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
        node_random_color = add_custom(
            group, 'MN_color_attribute_random', [-300, -150])

        link(node_input.outputs['Geometry'], node_color_set.inputs[0])
        link(node_random_color.outputs['Color'],
             node_color_common.inputs['Carbon'])
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
        node_animate.inputs['To Max'].default_value = len(
            coll_frames.objects) - 1

        link(to_animate.outputs[0], node_animate_frames.inputs[0])
        link(node_animate_frames.outputs[0], node_style.inputs[0])
        link(node_animate.outputs[0], node_animate_frames.inputs['Frame'])


def combine_join_geometry(group, node_list, output='Geometry', join_offset=300):
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
        link(named_att.outputs['Attribute'], node_split.inputs['Field'])
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

    transforms = utils.array_quaternions_from_dict(
        mol['biological_assemblies'])
    data_object = obj.create_data_object(
        array=transforms,
        name=f"data_assembly_{mol.name}"
    )
    tree_assembly = create_assembly_node_tree(
        name=mol.name,
        iter_list=mol['chain_ids'],
        data_object=data_object
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
    group = new_group(name=node_group_name)
    link = group.links.new

    n_assemblies = len(
        np.unique(obj.get_attribute(data_object, 'assembly_id')))

    node_group_instances = split_geometry_to_instances(
        name=f".MN_utils_split_{name}",
        iter_list=iter_list,
        attribute='chain_id'
    )

    node_group_assembly_instance = append('.MN_assembly_instance_chains')
    node_instances = add_custom(group, node_group_instances.name, [0, 0])
    node_assembly = add_custom(
        group, node_group_assembly_instance.name, [200, 0])
    node_assembly.inputs['data_object'].default_value = data_object

    out_sockets = outputs(group)
    out_sockets[list(out_sockets)[0]].name = "Instances"

    socket_info = (
        {"name": "Rotation",    "type": "NodeSocketFloat",
            "min": 0, "max": 1, "default": 1},
        {"name": "Translation", "type": "NodeSocketFloat",
            "min": 0, "max": 1, "default": 1},
        {"name": "assembly_id", "type": "NodeSocketInt",
            "min": 1, "max": n_assemblies, "default": 1}
    )

    for info in socket_info:
        socket = group.interface.items_tree.get(info['name'])
        if not socket:
            socket = group.interface.new_socket(
                info['name'], in_out='INPUT', socket_type=info['type'])
        socket.default_value = info['default']
        socket.min_value = info['min']
        socket.max_value = info['max']

        link(get_input(group).outputs[info['name']],
             node_assembly.inputs[info['name']])

    get_output(group).location = [400, 0]
    link(get_input(group).outputs[0], node_instances.inputs[0])
    link(node_instances.outputs[0], node_assembly.inputs[0])
    link(node_assembly.outputs[0], get_output(group).inputs[0])

    return group


def add_inverse_selection(group):
    output = get_output(group)
    if 'Inverted' not in output.inputs.keys():
        group.interface.new_socket(
            'Inverted', in_out='OUTPUT', socket_type='NodeSocketBool')

    loc = output.location
    bool_math = group.nodes.new("FunctionNodeBooleanMath")
    bool_math.location = [loc[0], -100]
    bool_math.operation = "NOT"

    group.links.new(
        output.inputs['Selection'].links[0].from_socket, bool_math.inputs[0])
    group.links.new(bool_math.outputs[0], output.inputs['Inverted'])


def custom_iswitch(
        name,
        iter_list,
        field='chain_id',
        dtype='BOOLEAN',
        default_values=None,
        prefix='',
        start=0
):
    """
    Creates a named `Index Switch` node. 

    Wraps an index switch node, giving the group names or each name in the `iter_list`.
    Uses the given field for the attribute name to use in the index switch, and optionally
    adds an offset value if the start value is non zero.

    If a list of default items is given, then it is recycled to fill the defaults for 
    each created socket in for the node.

    Parameters
    ----------
    name : str
        The name of the node group.
    iter_list : list
        The list of items to iterate over.
    field : str, optional
        The name of the attribute field. Defaults to 'chain_id'.
    default_values : list, optional
        The list of default values to assign to each item. Defaults to None.
    prefix : str, optional
        The prefix to add to the node names. Defaults to an empty string.
    start : int, optional
        The starting index for the node names. Defaults to 0.

    Returns
    -------
    group : bpy.types.NodeGroup
        The created node group.

    Raises
    ------
    NodeGroupCreationError
        If there was an error creating the node group.
    """
    iter_list = [str(i) for i in iter_list]
    group = bpy.data.node_groups.get(name)
    if group:
        return group

    socket_type = socket_types[dtype]
    group = new_group(name, geometry=False, fallback=False)

    # try creating the node group, otherwise on fail cleanup the created group and
    # report the error
    try:
        link = group.links.new
        node_input = get_input(group)
        node_output = get_output(group)
        node_attr = group.nodes.new('GeometryNodeInputNamedAttribute')
        node_attr.data_type = 'INT'
        node_attr.location = [0, 150]
        node_attr.inputs['Name'].default_value = str(field)

        node_iswitch = group.nodes.new('GeometryNodeIndexSwitch')
        node_iswitch.data_type = dtype

        link(node_attr.outputs['Attribute'], node_iswitch.inputs['Index'])

        # if there is as offset to the lookup values (say we want to start looking up
        # from 100 or 1000 etc) then we add a math node with that offset value
        if start != 0:
            node_math = group.nodes.new('ShaderNodeMath')
            node_math.operation = 'ADD'
            node_math.location = [0, 150]
            node_attr.location = [0, 300]

            node_math.inputs[1].default_value = start
            link(
                node_attr.outputs['Attribute'],
                node_math.inputs[0]
            )
            link(
                node_math.outputs['Value'],
                node_iswitch.inputs['Index']
            )

        # if there are custom values provided, create a dictionary lookup for those values
        # to assign to the sockets upon creation. If no default was given and the dtype
        # is colors, then generate a random pastel color for each value
        default_lookup = None
        if default_values:
            default_lookup = dict(zip(
                iter_list,
                itertools.cycle(default_values)
            ))
        elif dtype == 'RGBA':
            default_lookup = dict(zip(
                iter_list,
                [color.random_rgb() for i in iter_list]
            ))

        # for each item in the iter_list, we create a new socket on the interface for this
        # node group, and link it to the interface on the index switch. The index switch
        # currently starts with two items already, so once i > 1 we start to add
        # new items for the index switch as well
        for i, item in enumerate(iter_list):
            if i > 1:
                node_iswitch.index_switch_items.new()

            socket = group.interface.new_socket(
                name=f'{prefix}{item}',
                in_out='INPUT',
                socket_type=socket_type
            )
            #  if a set of default values was given, then use it for setting
            # the defaults on the created sockets of the node group
            if default_lookup:
                socket.default_value = default_lookup[item]

            link(
                node_input.outputs[socket.identifier],
                node_iswitch.inputs[str(i)]
            )

        socket_out = group.interface.new_socket(
            name='Color',
            in_out='OUTPUT',
            socket_type=socket_type
        )
        link(
            node_iswitch.outputs['Output'],
            node_output.inputs[socket_out.identifier]
        )

        return group

    # if something broke when creating the node group, delete whatever was created
    except Exception as e:
        node_name = group.name
        bpy.data.node_groups.remove(group)
        raise NodeGroupCreationError(
            f'Unable to make node group: {node_name}.\nError: {e}'
        )


def resid_multiple_selection(node_name, input_resid_string):
    """
    Returns a node group that takes an integer input and creates a boolean
    tick box for each item in the input list. Outputs are the selected
    residues and the inverse selection. Used for constructing chain
    selections in specific proteins.
    """

    # print(f'recieved input: {input_resid_string}')
    # do a cleanning of input string to allow fuzzy input from users
    for c in ";/+ .":
        if c in input_resid_string:
            input_resid_string = input_resid_string.replace(c, ',')

    for c in "_=:":
        if c in input_resid_string:
            input_resid_string = input_resid_string.replace(c, '-')

    # print(f'fixed input:{input_resid_string}')

    # parse input_resid_string into sub selecting string list
    sub_list = [item for item in input_resid_string.split(',') if item]

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

    prev = None
    for residue_id_index, residue_id in enumerate(sub_list):

        # add an new node of Select Res ID or MN_sek_res_id_range
        current_node = new_node("GeometryNodeGroup")

        # add an bool_math block
        bool_math = new_node("FunctionNodeBooleanMath")
        bool_math.location = [400, (residue_id_index+1) * node_sep_dis]
        bool_math.operation = "OR"

        if '-' in residue_id:
            # set two new inputs
            current_node.node_tree = append('MN_select_res_id_range')
            [resid_start, resid_end] = residue_id.split('-')[:2]
            socket_1 = residue_id_group.interface.new_socket(
                'res_id: Min', in_out='INPUT', socket_type='NodeSocketInt')
            socket_1.default_value = int(resid_start)
            socket_2 = residue_id_group.interface.new_socket(
                'res_id: Max', in_out='INPUT', socket_type='NodeSocketInt')
            socket_2.default_value = int(resid_end)

            # a residue range
            group_link(
                node_input.outputs[socket_1.identifier], current_node.inputs[0])
            group_link(
                node_input.outputs[socket_2.identifier], current_node.inputs[1])
        else:
            # create a node
            current_node.node_tree = append('MN_select_res_id_single')
            socket = residue_id_group.interface.new_socket(
                'res_id', in_out='INPUT', socket_type='NodeSocketInt')
            socket.default_value = int(residue_id)
            group_link(
                node_input.outputs[socket.identifier], current_node.inputs[0])

        # set the coordinates
        current_node.location = [200, (residue_id_index+1) * node_sep_dis]
        if not prev:
            # link the first residue selection to the first input of its OR block
            group_link(current_node.outputs['Selection'], bool_math.inputs[0])
        else:
            # if it is not the first residue selection, link the output to the previous or block
            group_link(current_node.outputs['Selection'], prev.inputs[1])

            # link the ouput of previous OR block to the current OR block
            group_link(prev.outputs[0], bool_math.inputs[0])
        prev = bool_math

    # add a output block
    residue_id_group_out = new_node("NodeGroupOutput")
    residue_id_group_out.location = [
        800, (residue_id_index + 1) / 2 * node_sep_dis]
    residue_id_group.interface.new_socket(
        'Selection', in_out='OUTPUT', socket_type='NodeSocketBool')
    residue_id_group.interface.new_socket(
        'Inverted', in_out='OUTPUT', socket_type='NodeSocketBool')
    group_link(prev.outputs[0], residue_id_group_out.inputs['Selection'])
    invert_bool_math = new_node("FunctionNodeBooleanMath")
    invert_bool_math.location = [
        600, (residue_id_index+1) / 3 * 2 * node_sep_dis]
    invert_bool_math.operation = "NOT"
    group_link(prev.outputs[0], invert_bool_math.inputs[0])
    group_link(invert_bool_math.outputs[0],
               residue_id_group_out.inputs['Inverted'])
    return residue_id_group
