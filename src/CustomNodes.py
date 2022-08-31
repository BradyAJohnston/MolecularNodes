import bpy

# creates the data block for a boolean chain node. Useful for when
# constructing custom selection nodes, to be able to repeatedly add the
# boolean selection node. If the data block already exists, it will return that 
# data block to be used, if it doesn't, it will create the required data block
# for a chained node.
def create_bool_chain_data():

    try:
        bool_chain_data = bpy.data.node_groups['MOL_utils_bool_chain']
    except:
        bool_chain_data = None
    
    # check to see if the data block already exists, if it does just return that instead.
    if bool_chain_data != None:
        return bool_chain_data
    else:
        bgd = bool_group_data = bpy.data.node_groups.new("MOL_utils_bool_chain", "GeometryNodeTree")



        bool_group_in = bgd.nodes.new("NodeGroupInput")
        bool_group_in.location = [-200, 0]

        bgd.inputs.new("NodeSocketInt", "number_chain_in")
        bgd.inputs.new("NodeSocketBool", "bool_chain_in")
        bgd.inputs.new("NodeSocketInt", "number_matched")
        bgd.inputs.new("NodeSocketBool", "bool_include")

        bool_group_out = bgd.nodes.new("NodeGroupOutput")
        bgd.outputs.new("NodeSocketInt", "number_chain_out")
        bgd.outputs.new("NodeSocketBool", "bool_chain_out")

        bool_group_out.location = [600, 0]

        node_num_comp = bgd.nodes.new("FunctionNodeCompare")
        node_num_comp.data_type = "INT"
        node_num_comp.operation = "EQUAL"
        node_num_comp.location = [100, 100]

        node_bool_math_1 = bgd.nodes.new("FunctionNodeBooleanMath")
        node_bool_math_1.location = [100, -100]

        node_bool_math_2 = bgd.nodes.new("FunctionNodeBooleanMath")
        node_bool_math_2.operation = "OR"
        node_bool_math_2.location = [400, 0]

        link = bool_group_data.links.new

        link(bool_group_in.outputs["number_chain_in"], node_num_comp.inputs[2])
        link(bool_group_in.outputs["number_chain_in"], bool_group_out.inputs["number_chain_out"])
        link(bool_group_in.outputs["bool_chain_in"], node_bool_math_2.inputs[0])
        link(bool_group_in.outputs["number_matched"], node_num_comp.inputs[3])
        link(bool_group_in.outputs["bool_include"], node_bool_math_1.inputs[1])
        link(node_num_comp.outputs[0], node_bool_math_1.inputs[0])
        link(node_bool_math_1.outputs[0], node_bool_math_2.inputs[1])
        link(node_bool_math_2.outputs[0], bool_group_out.inputs["bool_chain_out"])

        return bool_group_data


def add_bool_chain_node():
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


    try:
        bool_chain_data = bpy.data.node_groups['MOL_utils_bool_chain']
    except:
        bool_chain_data = None

    if bool_chain_data == None:
        bool_chain_data = create_bool_chain_data()
    

    # now that the data block is checked and setup, create a blank 
    # group node and link the data block to the node.
    new_node = node_mod.node_group.nodes.new("GeometryNodeGroup")
    new_node.node_tree = bool_chain_data

    return new_node

# add_bool_chain_node()

def create_node_group(node_name, input_list, label_prefix = "Chain "):
    """
    Given a an input_list, will create a node which takes an Integer input, 
    and has a boolean tick box for each item in the input list. The outputs will
    be the resulting selection and the inversion of the selection.
    Can contain a prefix for the resulting labels. Mostly used for constructing 
    chain selections when required for specific proteins.
    """

    # get the active object, might need to change to taking an object as an input
    # and making it active isntead, to be more readily applied to multiple objects

    obj = bpy.context.active_object


    # try to get the Molecular Nodes modifier and select it, if not create one and select it
    try:
        node_mod = obj.modifiers['MolecularNodes']
    except: 
        node_mod = None

    if node_mod == None:
        node_mod = obj.modifiers.new("MolecularNodes", "NODES")
        obj.modifiers.active = node_mod
    else:
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
    chain_number_node.inputs[0].default_value = 'chain_number'
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
        current_node.node_tree = create_bool_chain_data()
        current_node.location = [counter * node_sep_dis, 200]
        current_node.inputs["number_matched"].default_value = counter + 1
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
chain_list = chain_id_list
node_name = "MOL_" + str(output_name) + "_selection_chain"

# finally make the selection node!
create_node_group(
    node_name = node_name, 
    input_list = chain_list, 
    label_prefix = "Chain"
    )

# node = bpy.context.selected_nodes[0]
# bpy.ops.transform.translate('INVOKE_DEFAULT')