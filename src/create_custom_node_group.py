import bpy

# following this tutorial: https://www.youtube.com/watch?v=tQW2EnId1Ks

def create_node_group(context, operator, group_name, chain_names):

    # get the active object, might need to change to taking an object as an input
    # and making it active isntead, to be more readily applied to multiple objects

    obj = bpy.context.active_object


    # try to get the Molecular Nodes mofier, if not create one and select it
    # if already exists, just select it
    node_mod = obj.modifiers['MolecularNodes']
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

    bpy.ops.node.group_make()