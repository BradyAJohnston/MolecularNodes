import bpy

# try to get the Molecular Nodes modifier and select it, if not create one and select it
def create_starting_node_tree(collection_of_properties, obj):
    try:
        node_mod = obj.modifiers['MolecularNodes']
    except: 
        node_mod = None

    if node_mod == None:
        node_mod = obj.modifiers.new("MolecularNodes", "NODES")
        obj.modifiers.active = node_mod
    else:
        obj.modifiers.active = node_mod

    node_mod.node_group.name = "MOL_" + str(output_name)
    

    node_input = node_mod.node_group.nodes['Group Input']
    node_input.location = [-200, 0]
    node_output = node_mod.node_group.nodes['Group Output']
    node_output.location = [600, 0]

    # create an empty node group and link it to the atomic properties node group
    new_node_group = node_mod.node_group.nodes.new("GeometryNodeGroup")
    new_node_group.node_tree = bpy.data.node_groups["MOL_atomic_properties"]
    new_node_group.inputs['Properties'].default_value = collection_of_properties
    new_node_group.location = [0, 0]
    # resize the newly created node to be a bit wider
    node_mod.node_group.nodes[-1].width = 200



    colour_node_group = node_mod.node_group.nodes.new("GeometryNodeGroup")
    colour_node_group.node_tree = bpy.data.node_groups["MOL_style_colour"]
    colour_node_group.location = [300,0]
    node_mod.node_group.nodes[-1].width = 200

    link = node_mod.node_group.links.new

    link(node_input.outputs['Geometry'], new_node_group.inputs['Atoms'])
    link(new_node_group.outputs['Atoms'], colour_node_group.inputs['Atoms'])
    link(new_node_group.outputs['atomic_number'], colour_node_group.inputs['atomic_number'])
    link(colour_node_group.outputs['Atoms'], node_output.inputs['Geometry'])

    node_mod.node_group.outputs.new("NodeSocketColor", "Colour")
    link(colour_node_group.outputs['Colour'], node_output.inputs['Colour'])
    node_mod['Output_2_attribute_name'] = "Colour"

    mat = create_starting_material()

    colour_node_group.inputs['Material'].default_value = mat


def create_starting_material():
    try:
        mat = bpy.data.materials['MOL_atomic_material']
        return mat
    except:
        mat = None

    if mat == None:
        mat = bpy.data.materials.new('MOL_atomic_material')

    mat.use_nodes = True
    node_att = mat.node_tree.nodes.new("ShaderNodeAttribute")
    
    node_att.attribute_name = "Colour"
    node_att.location = [-300, 200]

    mat.node_tree.links.new(node_att.outputs['Color'], mat.node_tree.nodes['Principled BSDF'].inputs['Base Color'])

    return mat
    
# create_starting_material(base_model)
create_starting_node_tree(collection_of_properties=col_properties, obj = base_model)

base_model.select_set(True)
bpy.context.view_layer.objects.active = base_model

# bpy.context.active_object.data.materials.append(mat)