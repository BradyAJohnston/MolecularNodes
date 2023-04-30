import bpy
import numpy as np
from . import nodes

def get_transformations_pdbx(file_pdbx):
    import biotite.structure.io.pdbx as pdbx
    #The output transform_dict has an entry for each transformation, indexable by the string
    # integer of the assembly number (e.g. transform_dict.get('1')) which contains tuple of the 3x3 rotation 
    # matrix and the 1x3 transform matrix
    
    transform_dict = pdbx.convert._get_transformations(file_pdbx.get_category('pdbx_struct_oper_list'))
    
    return transform_dict

def get_transformations_pdb(file_pdb):
    from re import compile
    # get the lines where 'SMTRY' appears, which specify the actual symmetry operations
    sym_lines = np.array(file_pdb.lines)[np.char.rfind(np.array(file_pdb.lines), 'SMTRY') > 0]
    
    regex = compile('\d\.\d+') # find where there is a digit, a decimal, and then multiple more digits
    
    n_mat = int(len(sym_lines) / 4)
    sym_array = np.zeros([len(sym_lines),4], dtype=np.float32)
    
    
    for i in range(len(sym_lines)):
        sym_array[i] = np.array(regex.findall(sym_lines[i]), dtype = np.float32)
    
    transform_dict = []
    
    for i in range(n_mat):
        
        mat_start = i * 3
        mat_end = (i + 1) * 3
        
        transform_dict.append((
            sym_array[mat_start:mat_end, :3], 
            sym_array[mat_start:mat_end, 3:].reshape([1, 3])
        ))
    
    return transform_dict

def get_transformations_mmtf(all_assemblies, world_scale = 0.01):
    """
    returns a (N, 3, 4) numpy matrix, where N is the number of transformations required
    to build out the biological assembly. Currently only extracts and supports the 
    first biological assembly, but this should be straightforward to expand to more assemblies
    but would require more tweaking with the creation of the nodes themselves than this function
    """
    
    assembly = all_assemblies[0]['transformList']
    matrices = np.array([item['matrix'] for item in assembly]).reshape((len(assembly), 4, 4))
    return matrices

def create_assembly_node(name, trans_mat):
    
    node_mat = bpy.data.node_groups.get('MOL_RotTransMat_' + name)
    if node_mat:
        return node_mat
    
    node_mat = nodes.gn_new_group_empty('MOL_RotTransMat_' + name)
    node_mat.inputs.remove(node_mat.inputs['Geometry'])
    node_mat.nodes['Group Output'].location = [800, 0]
    node_mat.outputs['Geometry'].name = 'RotTransMat'
    
    node_transform_list = []
    for i, mat in enumerate(trans_mat):
        node = nodes.rotation_matrix(
            node_group=node_mat, 
            mat=mat, 
            location=[0, 0 - (300 * i)]
        )
        
        node_transform_list.append(node)
    
    node_transform_list.reverse()
    
    node_join = node_mat.nodes.new('GeometryNodeJoinGeometry')
    node_join.location = [300, 0]
    
    for node_transform in node_transform_list:
        node_mat.links.new(node_transform.outputs['Geometry'], node_join.inputs['Geometry'])
    
    node_mat.links.new(node_join.outputs['Geometry'], node_mat.nodes['Group Output'].inputs['RotTransMat'])
    
    return node_mat

def create_biological_assembly_node(name, transform_dict):
    
    node_bio = bpy.data.node_groups.get('MOL_assembly_' + name)
    if node_bio:
        return node_bio
    
    # try to create the assembly transformation nodes first, so 
    # if they fail, nothing else is created
    data_trans = create_assembly_node(name, transform_dict)
    
    node_bio = nodes.gn_new_group_empty('MOL_assembly_' + name)
    
    node_input = node_bio.nodes[bpy.app.translations.pgettext_data("Group Input",)]
    node_output = node_bio.nodes[bpy.app.translations.pgettext_data("Group Output",)]
    
    
    node_output.location = [400, 0]
    
    node_assembly = nodes.add_custom_node_group_to_node(node_bio, 'MOL_utils_bio_assembly', location=[0, 0])
    
    node_trans = nodes.add_custom_node_group_to_node(node_bio, data_trans.name, location = [-400, -200])
    
    link = node_bio.links.new
    
    link(node_input.outputs['Geometry'], node_assembly.inputs['Geometry'])
    link(node_trans.outputs['RotTransMat'], node_assembly.inputs['RotTransMat'])
    link(node_assembly.outputs['Instances'], node_output.inputs[0])
    
    inputs = (
        {'name': 'Scale Rotation', 
         'type': 'NodeSocketFloat', 
         'default': 1},
        {'name': 'Scale Translation', 
         'type': 'NodeSocketFloat', 
         'default': 1}
    )
    
    for input in inputs:
        name = input.get('name')
        type = input.get('type')
        default = input.get('default')
        
        node_bio.inputs.new(type, name)
        node_bio.inputs.get(name).default_value = default
        
        link(node_input.outputs[name], node_assembly.inputs[name])
    
    return node_bio
