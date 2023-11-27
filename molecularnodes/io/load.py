import bpy
import numpy as np
import warnings
import time

from ..blender import (
    coll, obj
)
from .. import data
from .. import color


bpy.types.Scene.MN_import_centre = bpy.props.BoolProperty(
    name = "Centre Structure", 
    description = "Move the imported Molecule on the World Origin",
    default = False
    )
bpy.types.Scene.MN_import_del_solvent = bpy.props.BoolProperty(
    name = "Remove Solvent", 
    description = "Delete the solvent from the structure on import",
    default = True
    )
bpy.types.Scene.MN_import_panel_selection = bpy.props.IntProperty(
    name = "MN_import_panel_selection", 
    description = "Import Panel Selection", 
    subtype = 'NONE',
    default = 0
)
bpy.types.Scene.MN_import_build_assembly = bpy.props.BoolProperty(
    name = 'Build Assembly', 
    default = False
)

class MolecularNodesObjectProperties(bpy.types.PropertyGroup):
    subframes: bpy.props.IntProperty(
        name = "Subframes", 
        description = "Number of subframes to interpolate for MD trajectories.", 
        default = 0
    )
    molecule_type: bpy.props.StringProperty(
        name = "Molecular Type", 
        description = "How the file was imported, dictating how MN interacts with it.", 
        default = ""
    )
    pdb_code: bpy.props.StringProperty(
        name = "PDB", 
        description = "PDB code used to download this structure.", 
        maxlen = 4,
        options={'HIDDEN'}
    )



def pdb_get_b_factors(file):
    """
    Get a list, which contains a numpy array for each model containing the b-factors.
    """
    b_factors = []
    for model in range(file.get_model_count()):
        atoms = file.get_structure(model = model + 1, extra_fields = ['b_factor'])
        b_factors.append(atoms.b_factor)
    return b_factors

def get_secondary_structure(array, file) -> np.array:
    """
    Gets the secondary structure annotation that is included in mmtf files and returns it as a numerical numpy array.

    Parameters:
    -----------
    array : numpy.array
        The molecular coordinates array, from mmtf.get_structure()
    file : mmtf.MMTFFile
        The MMTF file containing the secondary structure information, from mmtf.MMTFFile.read()

    Returns:
    --------
    atom_sse : numpy.array
        Numerical numpy array representing the secondary structure of the molecule.
    
    Description:
    ------------
    This function uses the biotite.structure package to extract the secondary structure information from the MMTF file.
    The resulting secondary structures are `1: Alpha Helix, 2: Beta-sheet, 3: loop`.
    """
    
    from biotite.structure import spread_residue_wise
    
    sec_struct_codes = {
        -1: "X",
        0 : "I",
        1 : "S",
        2 : "H",
        3 : "E",
        4 : "G",
        5 : "B",
        6 : "T",
        7 : "C"
    }
    
    dssp_to_abc = {
        "X" : 0,
        "I" : 1, #"a",
        "S" : 3, #"c",
        "H" : 1, #"a",
        "E" : 2, #"b",
        "G" : 1, #"a",
        "B" : 2, #"b",
        "T" : 3, #"c",
        "C" : 3 #"c"
    }
    
    try:
        sse = file["secStructList"]
    except KeyError:
        ss_int = np.full(len(array), 3)
        print('Warning: "secStructList" field missing from MMTF file. Defaulting \
            to "loop" for all residues.')
    else:
        ss_int = np.array(
            [dssp_to_abc.get(sec_struct_codes.get(ss)) for ss in sse], 
            dtype = int
        )
    atom_sse = spread_residue_wise(array, ss_int)
    
    return atom_sse


def comp_secondary_structure(array):
    """Use dihedrals to compute the secondary structure of proteins

    Through biotite built-in method derivated from P-SEA algorithm (Labesse 1997)
    Returns an array with secondary structure for each atoms where:
    - 0 = '' = non-protein or not assigned by biotite annotate_sse
    - 1 = a = alpha helix
    - 2 = b = beta sheet
    - 3 = c = coil

    Inspired from https://www.biotite-python.org/examples/gallery/structure/transketolase_sse.html
    """
    #TODO Port [PyDSSP](https://github.com/ShintaroMinami/PyDSSP)
    from biotite.structure import annotate_sse, spread_residue_wise

    conv_sse_char_int = {'a': 1, 'b': 2, 'c': 3, '': 0} 

    char_sse = annotate_sse(array)
    int_sse = np.array([conv_sse_char_int[char] for char in char_sse], dtype=int)
    atom_sse = spread_residue_wise(array, int_sse)
        
    return atom_sse

def create_molecule(array, 
                    name, 
                    centre = False, 
                    file = None,
                    calculate_ss = False,
                    del_solvent = False, 
                    style = 0,
                    collection = None, 
                    verbose = False
                    ):
    import biotite.structure as struc
    
    frames = None
    if isinstance(array, struc.AtomArrayStack):
        if array.stack_depth() > 1:
            frames = array
        array = array[0]
    
    # remove the solvent from the structure if requested
    if del_solvent:
        array = array[np.invert(struc.filter_solvent(array))]

    world_scale = 0.01
    locations = array.coord * world_scale
    
    centroid = np.array([0, 0, 0])
    if centre:
        centroid = struc.centroid(array) * world_scale
    

    # subtract the centroid from all of the positions to localise the molecule on the world origin
    if centre:
        locations = locations - centroid

    if not collection:
        collection = coll.mn()
    
    bonds_array = []
    bond_idx = []

    if array.bonds:
        bonds_array = array.bonds.as_array()
        bond_idx = bonds_array[:, [0, 1]]
        bond_types = bonds_array[:, 2].copy(order = 'C') # the .copy(order = 'C') is to fix a weird ordering issue with the resulting array
    
    mol = obj.create_object(name=name, collection=collection, locations=locations, edges=bond_idx)
    
    # Add information about the bond types to the model on the edge domain
    # Bond types: 'ANY' = 0, 'SINGLE' = 1, 'DOUBLE' = 2, 'TRIPLE' = 3, 'QUADRUPLE' = 4
    # 'AROMATIC_SINGLE' = 5, 'AROMATIC_DOUBLE' = 6, 'AROMATIC_TRIPLE' = 7
    # https://www.biotite-python.org/apidoc/biotite.structure.BondType.html#biotite.structure.BondType
    if array.bonds:
        obj.add_attribute(object = mol, name = 'bond_type', data = bond_types, type = "INT", domain = "EDGE")

    # The attributes for the model are initially defined as single-use functions. This allows
    # for a loop that attempts to add each attibute by calling the function. Only during this
    # loop will the call fail if the attribute isn't accessible, and the warning is reported
    # there rather than setting up a try: except: for each individual attribute which makes
    # some really messy code.
    
    # I still don't like this as an implementation, and welcome any cleaner approaches that 
    # anybody might have.
    
    def att_atomic_number():
        atomic_number = np.array(list(map(
            lambda x: data.elements.get(x, {'atomic_number': -1}).get("atomic_number"), 
            np.char.title(array.element))))
        return atomic_number
    
    def att_res_id():
        return array.res_id
    
    def att_res_name():
        other_res = []
        counter = 0
        id_counter = -1
        res_names = array.res_name
        res_ids = array.res_id
        res_nums  = []
        
        for name in res_names:
            res_num = data.residues.get(name, {'res_name_num': -1}).get('res_name_num')
            
            if res_num == 9999:
                if res_names[counter - 1] != name or res_ids[counter] != res_ids[counter - 1]:
                    id_counter += 1
                
                unique_res_name = str(id_counter + 100) + "_" + str(name)
                other_res.append(unique_res_name)
                
                num = np.where(np.isin(np.unique(other_res), unique_res_name))[0][0] + 100
                res_nums.append(num)
            else:
                res_nums.append(res_num)
            counter += 1

        mol['ligands'] = np.unique(other_res)
        return np.array(res_nums)

    
    def att_chain_id():
        chain_id = np.searchsorted(np.unique(array.chain_id), array.chain_id)
        return chain_id
    
    def att_entity_id():
        return array.entity_id
    
    def att_b_factor():
        return array.b_factor
    
    def att_vdw_radii():
        vdw_radii =  np.array(list(map(
            # divide by 100 to convert from picometres to angstroms which is what all of coordinates are in
            lambda x: data.elements.get(x, {'vdw_radii': 100}).get('vdw_radii', 100) / 100,  
            np.char.title(array.element)
            )))
        return vdw_radii * world_scale
    
    def att_atom_name():
        atom_name = np.array(list(map(
            lambda x: data.atom_names.get(x, -1), 
            array.atom_name
        )))
        
        return atom_name

    def att_lipophobicity():
        lipo = np.array(list(map(
            lambda x, y: data.lipophobicity.get(x, {"0": 0}).get(y, 0),
            array.res_name, array.atom_name
        )))
        
        return lipo
    
    def att_charge():
        charge = np.array(list(map(
            lambda x, y: data.atom_charge.get(x, {"0": 0}).get(y, 0),
            array.res_name, array.atom_name
        )))
        return charge
    
    def att_color():
        return color.color_chains(att_atomic_number(), att_chain_id()).reshape(-1)
    
    def att_is_alpha():
        return np.isin(array.atom_name, 'CA')
    
    def att_is_solvent():
        return struc.filter_solvent(array)
    
    def att_is_backbone():
        """
        Get the atoms that appear in peptide backbone or nucleic acid phosphate backbones.
        Filter differs from the Biotite's `struc.filter_peptide_backbone()` in that this
        includes the peptide backbone oxygen atom, which biotite excludes. Additionally 
        this selection also includes all of the atoms from the ribose in nucleic acids, 
        and the other phosphate oxygens.
        """
        backbone_atom_names = [
            'N', 'C', 'CA', 'O',                    # peptide backbone atoms
            "P", "O5'", "C5'", "C4'", "C3'", "O3'", # 'continuous' nucleic backbone atoms
            "O1P", "OP1", "O2P", "OP2",             # alternative names for phosphate O's
            "O4'", "C1'", "C2'", "O2'"              # remaining ribose atoms
        ]
        
        is_backbone = np.logical_and(
            np.isin(array.atom_name, backbone_atom_names), 
            np.logical_not(struc.filter_solvent(array))
        )
        return is_backbone
    
    def att_is_nucleic():
        return struc.filter_nucleotides(array)
    
    def att_is_peptide():
        aa = struc.filter_amino_acids(array)
        con_aa = struc.filter_canonical_amino_acids(array)
        
        return aa | con_aa
    
    def att_is_hetero():
        return array.hetero
    
    def att_is_carb():
        return struc.filter_carbohydrates(array)

    def att_sec_struct():
        if calculate_ss or not file:
            return comp_secondary_structure(array)
        else:
            return get_secondary_structure(array, file)
    

    

    
    # these are all of the attributes that will be added to the structure
    # TODO add capcity for selection of particular attributes to include / not include to potentially
    # boost performance, unsure if actually a good idea of not. Need to do some testing.
    attributes = (
        {'name': 'res_id',          'value': att_res_id,              'type': 'INT',     'domain': 'POINT'},
        {'name': 'res_name',        'value': att_res_name,            'type': 'INT',     'domain': 'POINT'},
        {'name': 'atomic_number',   'value': att_atomic_number,       'type': 'INT',     'domain': 'POINT'},
        {'name': 'b_factor',        'value': att_b_factor,            'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'vdw_radii',       'value': att_vdw_radii,           'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'chain_id',        'value': att_chain_id,            'type': 'INT',     'domain': 'POINT'},
        {'name': 'entity_id',       'value': att_entity_id,           'type': 'INT',     'domain': 'POINT'},
        {'name': 'atom_name',       'value': att_atom_name,           'type': 'INT',     'domain': 'POINT'},
        {'name': 'lipophobicity',   'value': att_lipophobicity,       'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'charge',          'value': att_charge,              'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'Color',           'value': att_color,               'type': 'FLOAT_COLOR',   'domain': 'POINT'},
        
        {'name': 'is_backbone',     'value': att_is_backbone,         'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_alpha_carbon', 'value': att_is_alpha,            'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_solvent',      'value': att_is_solvent,          'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_nucleic',      'value': att_is_nucleic,          'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_peptide',      'value': att_is_peptide,          'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_hetero',       'value': att_is_hetero,           'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_carb',         'value': att_is_carb,             'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'sec_struct',      'value': att_sec_struct,          'type': 'INT',     'domain': 'POINT'}
    )
    
    # assign the attributes to the object
    for att in attributes:
        if verbose:
            start = time.process_time()
        try:
            obj.add_attribute(mol, att['name'], att['value'](), att['type'], att['domain'])
            if verbose:
                print(f'Added {att["name"]} after {time.process_time() - start} s')
        except :
            if verbose:
                warnings.warn(f"Unable to add attribute: {att['name']}")
                print(f'Failed adding {att["name"]} after {time.process_time() - start} s')

    if frames:
        try:
            b_factors = pdb_get_b_factors(file)
        except:
            b_factors = None
        
        coll_frames = coll.frames(mol.name, parent = coll.data())
        
        for i, frame in enumerate(frames):
            obj_frame = obj.create_object(
                name = mol.name + '_frame_' + str(i), 
                collection=coll_frames, 
                locations= frame.coord * world_scale - centroid
            )
            if b_factors:
                try:
                    obj.add_attribute(obj_frame, 'b_factor', b_factors[i])
                except:
                    b_factors = False
        
        # disable the frames collection so it is not seen
        # bpy.context.view_layer.layer_collection.children[''].children[coll_frames.name].exclude = True
    else:
        coll_frames = None
    
    mol.mn['molcule_type'] = 'pdb'
    
    # add custom properties to the actual blender object, such as number of chains, biological assemblies etc
    # currently biological assemblies can be problematic to holding off on doing that
    try:
        mol['chain_id_unique'] = list(np.unique(array.chain_id))
    except:
        warnings.warn('No chain information detected.')
    
    try: 
        mol['entity_names'] = [ent['description'] for ent in file['entityList']]
    except:
        pass
    
    return mol, coll_frames



