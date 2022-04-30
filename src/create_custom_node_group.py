import bpy

# following this tutorial: https://www.youtube.com/watch?v=tQW2EnId1Ks

def create_node_group(node_name):

    # get the active object, might need to change to taking an object as an input
    # and making it active isntead, to be more readily applied to multiple objects

    obj = bpy.context.active_object


    # try to get the Molecular Nodes mofier, if not create one and select it
    # if already exists, just select it
    try:
        node_mod = obj.modifiers['MolecularNodes']
    except: 
        node_mod = None

    if node_mod == None:
        node_mod = obj.modifiers.new("MolecularNodes", "NODES")
        obj.modifiers.active = node_mod
    else:
        obj.modifiers.active = node_mod

    link = node_mod.node_group.links.new

    node_output = node_mod.node_group.nodes['Group Output']
    node_input = node_mod.node_group.nodes['Group Input']
    node_output.location = (800, 0)

    # create the node that will join all of the geometry together
    join_geometry = node_mod.node_group.nodes.new("GeometryNodeJoinGeometry")
    join_geometry.location = (500, 0)

    link(node_input.outputs['Geometry'], join_geometry.inputs['Geometry'])
    link(join_geometry.outputs['Geometry'], node_output.inputs['Geometry'])

    # chain_group = node_mod.node_group.nodes.new("GeometryNodeGroup")
    
    chain_group = bpy.data.node_groups.new(node_name, "GeometryNodeTree")
    chain_group_in = chain_group.nodes.new("NodeGroupInput")
    chain_group_in.location = [-200, 0]
    chain_group.inputs.new("NodeSocketBool", "Chain A")
    chain_group.inputs.new("NodeSocketBool", "Chain B")
    new_node = chain_group.nodes.new

    chain_group_out = chain_group.nodes.new("NodeGroupOutput")
    chain_group_out.location = [200, 0]

    chain_group.outputs.new("NodeSocketBool", "Selection")

    node_join = new_node("GeometryNodeJoinGeometry")

    # node_input = chain_group.nodes['Group Input']
    # node_output = chain_group.nodes['Group Output']

    # link = chain_group.node_group.nodes.link
    # link(node_input.outputs[0], node_join.inputs[0])
    # link(node_join.outputs[0], node_output.inputs[0])

    new_node_group = node_mod.node_group.nodes.new("GeometryNodeGroup")

    new_node_group.node_tree = bpy.data.node_groups[chain_group.name]

create_node_group("New Test Node group")