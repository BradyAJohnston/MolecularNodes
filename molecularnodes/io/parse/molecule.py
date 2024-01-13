from abc import ABCMeta
import warnings
import time
import numpy as np
import bpy

from ...blender import nodes, coll, obj
from ... import utils, data, color


class Molecule(metaclass=ABCMeta):
    """
    Abstract base class for representing a molecule.

    Methods
    -------
    entity_ids(as_int=False)
        Get the unique entity IDs of the molecule.
    chain_ids(as_int=False)
        Get the unique chain IDs of the molecule.
    create_model(name='NewMolecule', style='spheres', build_assembly=False, centre=False, del_solvent=True, collection=None, verbose=False)
        Create a model in the 3D scene.

    Attributes
    ----------
    structure : ndarray
        The structure of the molecule.
    """

    def chain_ids(self, as_int=False):
        """
        Get the unique chain IDs of the molecule.

        Parameters
        ----------
        as_int : bool, optional
            Whether to return the chain IDs as integers. Default is False.

        Returns
        -------
        ndarray
            The unique chain IDs of the molecule.
        """
        if as_int:
            return np.unique(self.structure.chain_id, return_inverse=True)[1]
        return np.unique(self.structure.chain_id)

    def create_model(
        self,
        name: str = 'NewMolecule',
        style: str = 'spheres',
        selection: np.ndarray = None,
        build_assembly=False,
        centre: bool = False,
        del_solvent: bool = True,
        collection=None,
        verbose: bool = True,
    ) -> bpy.types.Object:
        """
        Create a model in the 3D scene.

        Parameters
        ----------
        name : str, optional
            The name of the model. Default is 'NewMolecule'.
        centre : bool, optional
            Whether to center the model in the scene. Default is False.
        style : str, optional
            The style of the model. Default is 'spheres'.
        del_solvent : bool, optional
            Whether to delete solvent molecules. Default is True.
        collection : str, optional
            The collection to add the model to. Default is None.
        verbose : bool, optional
            Whether to print verbose output. Default is False.

        Returns
        -------
        model : bpy.types.Object
            The created 3D model, as an object in the 3D scene.
        """
        from biotite import InvalidFileError

        if selection:
            array = self.structure[selection]
        else:
            array = self.structure

        model, frames = _create_model(
            array=array,
            name=name,
            centre=centre,
            del_solvent=del_solvent,
            style=style,
            collection=collection,
            verbose=verbose,
        )

        if style:
            nodes.create_starting_node_tree(
                object=model,
                coll_frames=frames,
                style=style,
            )

        try:
            model['entity_ids'] = self.entity_ids
        except AttributeError:
            model['entity_ids'] = None

        try:
            model['biological_assemblies'] = self.assemblies()
        except InvalidFileError:
            pass

        if build_assembly and style:
            nodes.assembly_insert(model)

        return model

    def assemblies(self, as_array=False):
        from biotite import InvalidFileError
        try:
            assemblies_info = self._assemblies()
        except InvalidFileError:
            return None

        if isinstance(assemblies_info, dict) and as_array:
            return utils.array_quaternions_from_dict(assemblies_info)

        return assemblies_info


def _create_model(array,
                  name=None,
                  centre=False,
                  del_solvent=False,
                  style='spherers',
                  collection=None,
                  world_scale=0.01,
                  verbose=True
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
        # the .copy(order = 'C') is to fix a weird ordering issue with the resulting array
        bond_types = bonds_array[:, 2].copy(order='C')

    mol = obj.create_object(name=name, collection=collection,
                            vertices=locations, edges=bond_idx)

    # Add information about the bond types to the model on the edge domain
    # Bond types: 'ANY' = 0, 'SINGLE' = 1, 'DOUBLE' = 2, 'TRIPLE' = 3, 'QUADRUPLE' = 4
    # 'AROMATIC_SINGLE' = 5, 'AROMATIC_DOUBLE' = 6, 'AROMATIC_TRIPLE' = 7
    # https://www.biotite-python.org/apidoc/biotite.structure.BondType.html#biotite.structure.BondType
    if array.bonds:
        obj.set_attribute(object=mol, name='bond_type',
                          data=bond_types, type="INT", domain="EDGE")

    # The attributes for the model are initially defined as single-use functions. This allows
    # for a loop that attempts to add each attibute by calling the function. Only during this
    # loop will the call fail if the attribute isn't accessible, and the warning is reported
    # there rather than setting up a try: except: for each individual attribute which makes
    # some really messy code.

    # I still don't like this as an implementation, and welcome any cleaner approaches that
    # anybody might have.

    def att_atomic_number():
        atomic_number = np.array([
            data.elements.get(x, {'atomic_number': -1}).get("atomic_number")
            for x in np.char.title(array.element)
        ])
        return atomic_number

    def att_atom_id():
        return array.atom_id

    def att_res_id():
        return array.res_id

    def att_res_name():
        other_res = []
        counter = 0
        id_counter = -1
        res_names = array.res_name
        res_ids = array.res_id
        res_nums = []

        for name in res_names:
            res_num = data.residues.get(
                name, {'res_name_num': -1}).get('res_name_num')

            if res_num == 9999:
                if res_names[counter - 1] != name or res_ids[counter] != res_ids[counter - 1]:
                    id_counter += 1

                unique_res_name = str(id_counter + 100) + "_" + str(name)
                other_res.append(unique_res_name)

                num = np.where(np.isin(np.unique(other_res), unique_res_name))[
                    0][0] + 100
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

    def att_occupancy():
        return array.occupancy

    def att_vdw_radii():
        vdw_radii = np.array(list(map(
            # divide by 100 to convert from picometres to angstroms which is what all of coordinates are in
            lambda x: data.elements.get(
                x, {'vdw_radii': 100}).get('vdw_radii', 100) / 100,
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
        return color.color_chains(att_atomic_number(), att_chain_id())

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
            "P", "O5'", "C5'", "C4'", "C3'", "O3'",  # 'continuous' nucleic backbone atoms
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
        return array.sec_struct

    # these are all of the attributes that will be added to the structure
    # TODO add capcity for selection of particular attributes to include / not include to potentially
    # boost performance, unsure if actually a good idea of not. Need to do some testing.
    attributes = (
        {'name': 'res_id',          'value': att_res_id,
            'type': 'INT',     'domain': 'POINT'},
        {'name': 'res_name',        'value': att_res_name,
            'type': 'INT',     'domain': 'POINT'},
        {'name': 'atomic_number',   'value': att_atomic_number,
            'type': 'INT',     'domain': 'POINT'},
        {'name': 'b_factor',        'value': att_b_factor,
            'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'occupancy',       'value': att_occupancy,
            'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'vdw_radii',       'value': att_vdw_radii,
            'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'chain_id',        'value': att_chain_id,
            'type': 'INT',     'domain': 'POINT'},
        {'name': 'entity_id',       'value': att_entity_id,
            'type': 'INT',     'domain': 'POINT'},
        {'name': 'atom_id',         'value': att_atom_id,
            'type': 'INT',     'domain': 'POINT'},
        {'name': 'atom_name',       'value': att_atom_name,
            'type': 'INT',     'domain': 'POINT'},
        {'name': 'lipophobicity',   'value': att_lipophobicity,
            'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'charge',          'value': att_charge,
            'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'Color',           'value': att_color,
            'type': 'FLOAT_COLOR',   'domain': 'POINT'},

        {'name': 'is_backbone',     'value': att_is_backbone,
            'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_alpha_carbon', 'value': att_is_alpha,
            'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_solvent',      'value': att_is_solvent,
            'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_nucleic',      'value': att_is_nucleic,
            'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_peptide',      'value': att_is_peptide,
            'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_hetero',       'value': att_is_hetero,
            'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_carb',         'value': att_is_carb,
            'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'sec_struct',      'value': att_sec_struct,
            'type': 'INT',     'domain': 'POINT'}
    )

    # assign the attributes to the object
    for att in attributes:
        if verbose:
            start = time.process_time()
        # try:
        obj.set_attribute(
            mol, att['name'], att['value'](), att['type'], att['domain'])
        if verbose:
            print(
                f'Added {att["name"]} after {time.process_time() - start} s')
        # except:
        #     if verbose:
        #         warnings.warn(f"Unable to add attribute: {att['name']}")
        #         print(
        #             f'Failed adding {att["name"]} after {time.process_time() - start} s')

    coll_frames = None
    if frames:
        coll_frames = coll.frames(mol.name, parent=coll.data())
        for i, frame in enumerate(frames):
            obj.create_object(
                name=mol.name + '_frame_' + str(i),
                collection=coll_frames,
                vertices=frame.coord * world_scale - centroid
            )

    mol.mn['molcule_type'] = 'pdb'

    # add custom properties to the actual blender object, such as number of chains, biological assemblies etc
    # currently biological assemblies can be problematic to holding off on doing that
    try:
        mol['chain_ids'] = list(np.unique(array.chain_id))
    except AttributeError:
        mol['chain_ids'] = None
        warnings.warn('No chain information detected.')

    return mol, coll_frames
