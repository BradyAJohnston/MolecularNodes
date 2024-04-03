from typing import Optional, Any, Union, List, Tuple, Dict
import warnings
import time
import numpy as np
from functools import singledispatchmethod

import biotite.structure as struc
from biotite import InvalidFileError

import bpy
from bpy.app.handlers import persistent

from . import coll, obj, nodes

from ..io.parse.molecule import MoleculeAtomArray, AtomList
from ..io.parse.mda import MDA

from .obj import set_attribute, get_attribute, create_object
from .nodes import assembly_insert, create_starting_node_tree
from .obj import evaluated as evaluate, create_object
from .coll import mn as coll_mn, data as coll_data, frames as coll_frames

from ..data import elements, residues, atom_names, lipophobicity, atom_charge
from ..color import color_chains

from ..utils import lerp

from ..pkg import start_logging

class MoleculeInBlender:
    """
    Represents a molecule in Blender.

    Converts a given molecule into a Blender object.
    Long term goal to use molecule as data twin. So that any changes in molecule will be reflected in blender object.
    
    Attributes:
    ----
    name (str): Name of the molecule.
    molecule (bpy.types.Object): 3d object in blender. 
    frame_collection (bpy.types.Collection) : s. issue #454
    """

    __slots__ = ["name", "object", "frames", "universe_reps", "atom_reps", "in_memory", "log"]

    def __init__(self, name : str, object : bpy.types.Object, frames : bpy.types.Collection, log, in_memory : bool = False, ):

        if not log:
            log = start_logging(logfile_name=name)

        # if the session already exists, load the existing session
        if hasattr(bpy.types.Scene, "mda_session"):
            warnings.warn("The existing mda session is loaded.")
            log.warning("The existing mda session is loaded.")
            existing_session = bpy.types.Scene.mda_session
            self.__dict__ = existing_session.__dict__
            return

        self.name:str = name
        self.object: bpy.types.Object = object #TODO: rename object, so its function is clearer
        self.frames: bpy.types.Collection = frames #TODO: rename frames, so its function is clearer s. Issue #454
        self.universe_reps : Dict = dict()
        self.atom_reps : Dict = dict()
        self.in_memory = in_memory
        self.log = log

        if self.in_memory:
            return
        bpy.types.Scene.mda_session = self
        bpy.app.handlers.frame_change_post.append(
            self._update_trajectory_handler_wrapper()
        )
        bpy.app.handlers.depsgraph_update_pre.append(
            self._update_style_handler_wrapper()
        )
        self.log.info("MDAnalysis session is initialized.")


    @singledispatchmethod
    @classmethod
    def from_molecule(cls, molecule, **kwargs) -> "MoleculeInBlender":
        raise NotImplementedError(f"Class {molecule.__class__=} is not registrated")
    
    @from_molecule.register
    @classmethod
    def _(cls, molecule: MoleculeAtomArray, 
                      style: str = 'spheres', 
                      selection: np.ndarray = None, 
                      build_assembly: bool = False, 
                      centre : bool = False, 
                      del_solvent: bool= True, 
                      collection : bpy.types.Collection = None, 
                      verbose: bool = False
                    ):
        
        if selection:
            array = molecule[selection]._atoms
        else:
            array = molecule._atoms

        model, frames = _create_model(
            array=array,
            name=molecule.name,
            centre=centre,
            del_solvent=del_solvent,
            style=style,
            collection=collection,
            verbose=verbose,
        )

        if style:
            create_starting_node_tree(object=model,
                                      coll_frames=frames,
                                      style=style,
                                    )
            
        try:
            model['entity_ids'] = molecule.entity_ids
        except AttributeError:
            model['entity_ids'] = None

        try:
            model['biological_assemblies'] = molecule.assemblies()
        except InvalidFileError:
            pass

        if build_assembly and style:
            assembly_insert(model)

        return cls(name=molecule.name, object=model, frames=frames)
    
    @from_molecule.register
    @classmethod
    def _(cls, molecule:MDA, 
        style: str = "vdw", 
        selection: str = "all",
        name: str = "atoms",
        custom_selections: Dict[str, str] = {},
        frame_mapping: np.ndarray = None,
        subframes: int = 0,
        in_memory: bool = False
    ):
        """

        mix out of init and show
        """
        log = start_logging(logfile_name=molecule.name)

        #if in_memory:
        #    mol_object = cls.in_memory(
        #        atoms=molecule._atoms,
        #        style=style,
        #        selection=selection,
        #        name=molecule.name,
        #        custom_selections=custom_selections
        #    )
        #    if frame_mapping is not None:
        #        warnings.warn("Custom frame_mapping not supported"
        #                      "when in_memory is on.")
        #    if subframes != 0:
        #        warnings.warn("Custom subframes not supported"
        #                      "when in_memory is on.")
        #    log.info(f"{atoms} is loaded in memory.")
        #    return cls(name=molecule.name, object = mol_object, frames=None, log=log)
        #
        #if isinstance(atoms, mda.Universe):
        #    atoms = atoms.select_atoms(selection)

        universe = molecule.universe

        # if any frame_mapping is out of range, then raise an error
        if frame_mapping and (len(frame_mapping) > universe.trajectory.n_frames):
            raise ValueError("one or more mapping values are"
                             "out of range for the trajectory")
        

        mol_object =_process_atomgroup(
            ag=molecule,
            frame_mapping=frame_mapping,
            subframes=subframes,
            style=style,
            return_object=True)

        # add the custom selections if they exist
        #for sel_name, sel in custom_selections.items():
        #    try:
        #        ag = universe.select_atoms(sel)
        #        if ag.n_atoms == 0:
        #            raise ValueError("Selection is empty")
        #        self._process_atomgroup(
        #            ag=molecule._atoms,
        #            frame_mapping=frame_mapping,
        #            subframes=subframes,
        #            name=sel_name,
        #            style=style,
        #            return_object=False
        #        )
        #    except ValueError:
        #        warnings.warn(
        #            "Unable to add custom selection: {}".format(name))

        bpy.context.view_layer.objects.active = mol_object
        log.info(f"{len(molecule)} is loaded.")

        return cls(name=molecule.name, object = mol_object, frames=None, log=log)

    def __repr__(self) -> str:
        return f"<MoleculeInBlender: {self.name}>"
    
    def __setattr__(self, name: str, value: Any) -> None:
        """
        Set attr to the instance of the class, if the attr is in the __slots__, 
        else to the molecule object of the instance.

        Warning: Does not work as it should.
        """

        if name in self.__slots__:
            object.__setattr__(self, name, value )
        else:
            set_attribute(object=self.object, name=name, data=value)
    
    def __getattr__(self, name):
        """
        Get attr from the instance's obejct of the class, if the attr was not found in the __slots__ by the __getattribute__ method.

        Raises:
        ---
        AttributeError: If the attribute is not found in the instance of the class nor in the molecule object.
        """

        return get_attribute(object=self.object, name=name)
    
    def list_attributes(self, evaluate=False) -> list | None:
        """
        Returns a list of attribute names for the object.

        Parameters
        ----------
        evaluate : bool, optional
            Whether to first evaluate the modifiers on the object before listing the 
            available attributes.

        Returns
        -------
        list[str] | None
            A list of attribute names if the molecule object exists, None otherwise.
        """
        if not self.object:
            warnings.warn("No object created")
            return None
        if evaluate:
            return list(evaluate(self.object).data.attributes.keys())

        return list(self.object.data.attributes.keys())
    
    @persistent
    def _update_trajectory(self, frame):
        """
        The function that will be called when the frame changes.
        It will update the positions and selections of the atoms in the scene.
        """
        for rep_name in self.rep_names:
            universe = self.universe_reps[rep_name]["universe"]
            frame_mapping = self.universe_reps[rep_name]["frame_mapping"]
            subframes = bpy.data.objects[rep_name].mn['subframes']

            if frame < 0:
                continue

            if frame_mapping:
                # add the subframes to the frame mapping
                frame_map = np.repeat(frame_mapping, subframes + 1)
                # get the current and next frames
                frame_a = frame_map[frame]
                frame_b = frame_map[frame + 1]

            else:
                # get the initial frame
                if subframes == 0:
                    frame_a = frame
                else:
                    frame_a = int(frame / (subframes + 1))

                # get the next frame
                frame_b = frame_a + 1

            if frame_a >= universe.trajectory.n_frames:
                continue

            ag_rep = self.atom_reps[rep_name]
            mol_object = bpy.data.objects[rep_name]

            # set the trajectory at frame_a
            universe.trajectory[frame_a]

            if subframes > 0:
                fraction = frame % (subframes + 1) / (subframes + 1)

                # get the positions for the next frame
                locations_a = ag_rep.positions

                if frame_b < universe.trajectory.n_frames:
                    universe.trajectory[frame_b]
                locations_b = ag_rep.positions

                # interpolate between the two sets of positions
                locations = lerp(locations_a, locations_b, t=fraction)
            else:
                locations = ag_rep.positions

            # if the class of AtomGroup is UpdatingAtomGroup
            # then update as a new mol_object
            if isinstance(ag_rep.ag, mda.core.groups.UpdatingAtomGroup):
                mol_object.data.clear_geometry()
                mol_object.data.from_pydata(
                    ag_rep.positions,
                    ag_rep.bonds,
                    faces=[])
                for att_name, att in ag_rep._attributes_2_blender.items():
                    obj.set_attribute(
                        mol_object, att_name, att["value"], att["type"], att["domain"]
                    )
                mol_object['chain_id'] = ag_rep.chain_ids
                mol_object['atom_type_unique'] = ag_rep.atom_type_unique
                mol_object.mn['subframes'] = subframes
            else:
                # update the positions of the underlying vertices
                obj.set_attribute(mol_object, 'position', locations)

    @persistent
    def _update_trajectory_handler_wrapper(self):
        """
        A wrapper for the update_trajectory function because Blender
        requires the function to be taking one argument.
        """
        def update_trajectory_handler(scene):
            frame = scene.frame_current
            self._update_trajectory(frame)

        return update_trajectory_handler

    @persistent
    def _update_style_handler_wrapper(self):
        """
        A wrapper for the update_style function because Blender
        requires the function to be taking one argument.
        """
        def update_style_handler(scene):
            self._remove_deleted_mol_objects()
            # TODO: check for topology changes
            # TODO: update for style changes

        return update_style_handler

    @persistent
    def _remove_deleted_mol_objects(self):
        """
        Remove the deleted mol objects (e.g. due to operations inside Blender)
        from the session.
        """
        for rep_name in self.rep_names:
            if rep_name not in bpy.data.objects:
                self.rep_names.remove(rep_name)
                del self.atom_reps[rep_name]
                del self.universe_reps[rep_name]

    def _dump(self, blender_save_loc):
        """
        Dump the session as a pickle file 
        """
        log = start_logging(logfile_name="mda")
        # get blender_save_loc
        blender_save_loc = blender_save_loc.split(".blend")[0]
        with open(f"{blender_save_loc}.mda_session", "wb") as f:
            pickle.dump(self, f)
        log.info("MDAnalysis session is dumped to {}".
                    format(blender_save_loc))

    @classmethod
    def _rejuvenate(cls, mol_objects):
        """
        Rejuvenate the session from a pickle file in the default location
        (`~/.blender_mda_session/`).
        """
        log = start_logging(logfile_name="mda")

        # get session name from mol_objects dictionary
        blend_file_name = bpy.data.filepath.split(".blend")[0]
        try:
            with open(f"{blend_file_name}.mda_session", "rb") as f:
                cls = pickle.load(f)
        except FileNotFoundError:
            return None
        bpy.app.handlers.frame_change_post.append(
            cls._update_trajectory_handler_wrapper()
        )
        bpy.app.handlers.depsgraph_update_pre.append(
            cls._update_style_handler_wrapper()
        )
        log.info("MDAnalysis session is loaded from {}".
                    format(blend_file_name))
        return cls

    

    

def _create_model(array,
                  name=None,
                  centre=False,
                  del_solvent=False,
                  style='spherers',
                  collection=None,
                  world_scale=0.01,
                  verbose=False
                  ) -> Tuple[bpy.types.Object, bpy.types.Collection]:
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
        collection = coll_mn()


    bonds_array = []
    bond_idx = []

    if array.bonds:
        bonds_array = array.bonds.as_array()
        bond_idx = bonds_array[:, [0, 1]]
        # the .copy(order = 'C') is to fix a weird ordering issue with the resulting array
        bond_types = bonds_array[:, 2].copy(order='C')

    mol = create_object(name=name, collection=collection,
                               vertices=locations, edges=bond_idx)

    # Add information about the bond types to the model on the edge domain
    # Bond types: 'ANY' = 0, 'SINGLE' = 1, 'DOUBLE' = 2, 'TRIPLE' = 3, 'QUADRUPLE' = 4
    # 'AROMATIC_SINGLE' = 5, 'AROMATIC_DOUBLE' = 6, 'AROMATIC_TRIPLE' = 7
    # https://www.biotite-python.org/apidoc/biotite.structure.BondType.html#biotite.structure.BondType
    if array.bonds:
        set_attribute(mol, name='bond_type', data=bond_types,
                             type="INT", domain="EDGE")

    # The attributes for the model are initially defined as single-use functions. This allows
    # for a loop that attempts to add each attibute by calling the function. Only during this
    # loop will the call fail if the attribute isn't accessible, and the warning is reported
    # there rather than setting up a try: except: for each individual attribute which makes
    # some really messy code.

    # I still don't like this as an implementation, and welcome any cleaner approaches that
    # anybody might have.

    def att_atomic_number():
        atomic_number = np.array([
            elements.get(x, {'atomic_number': -1}).get("atomic_number")
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
            res_num = residues.get(
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
        return np.unique(array.chain_id, return_inverse=True)[1]

    def att_entity_id():
        return array.entity_id

    def att_b_factor():
        return array.b_factor

    def att_occupancy():
        return array.occupancy

    def att_vdw_radii():
        vdw_radii = np.array(list(map(
            # divide by 100 to convert from picometres to angstroms which is what all of coordinates are in
            lambda x: elements.get(
                x, {'vdw_radii': 100}).get('vdw_radii', 100) / 100,
            np.char.title(array.element)
        )))
        return vdw_radii * world_scale

    def att_atom_name():
        atom_name = np.array(list(map(
            lambda x: atom_names.get(x, -1),
            array.atom_name
        )))

        return atom_name

    def att_lipophobicity():
        lipo = np.array(list(map(
            lambda x, y: lipophobicity.get(x, {"0": 0}).get(y, 0),
            array.res_name, array.atom_name
        )))

        return lipo

    def att_charge():
        charge = np.array(list(map(
            lambda x, y: atom_charge.get(x, {"0": 0}).get(y, 0),
            array.res_name, array.atom_name
        )))
        return charge

    def att_color():
        return color_chains(att_atomic_number(), att_chain_id())

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
        try:
            set_attribute(mol, name=att['name'], data=att['value'](
            ), type=att['type'], domain=att['domain'])
            if verbose:
                print(
                    f'Added {att["name"]} after {time.process_time() - start} s')
        except:
            if verbose:
                warnings.warn(f"Unable to add attribute: {att['name']}")
                print(
                    f'Failed adding {att["name"]} after {time.process_time() - start} s')

    coll_frames = None
    if frames:

        coll_frames = coll_frames(mol.name, parent=coll_data())
        for i, frame in enumerate(frames):
            frame = create_object(
                name=mol.name + '_frame_' + str(i),
                collection=coll_frames,
                vertices=frame.coord * world_scale - centroid
            )
            # TODO if update_attribute
            # bl.obj.set_attribute(attribute)

    mol.mn['molcule_type'] = 'pdb'

    # add custom properties to the actual blender object, such as number of chains, biological assemblies etc
    # currently biological assemblies can be problematic to holding off on doing that
    try:
        mol['chain_ids'] = list(np.unique(array.chain_id))
    except AttributeError:
        mol['chain_ids'] = None
        warnings.warn('No chain information detected.')

    return mol, coll_frames

def _process_atomgroup(
        ag : MDA,
        world_scale=0.01,
        frame_mapping=None,
        subframes=0,
        style="vdw",
        node_setup=True,
        return_object=False,
    ):
        """
        process the atomgroup in the Blender scene.

        Parameters:
        ----------
        ag : MDAnalysis.AtomGroup
            The atomgroup to add in the scene.
        frame_mapping : np.ndarray
            The frame mapping for the trajectory in Blender frame indices. Default: None
        subframes : int
            The number of subframes to interpolate between each frame.
        name : str
            The name of the atomgroup. Default: 'atoms'
        style : str
            The style of the atoms. Default: 'vdw'
        node_setup : bool
            Whether to add the node tree for the atomgroup. Default: True
        return_object : bool
            Whether to return the blender object or not. Default: False
        """
        # create the initial model
        mol_object = obj.create_object(
            name=ag.name,
            collection=coll.mn(),
            vertices=ag.positions,
            edges=ag.bonds,
        )

        # add the attributes for the model in blender
        for att_name, att in ag._attributes_2_blender.items():
            obj.set_attribute(
                mol_object, att_name, att["value"], att["type"], att["domain"]
            )
        mol_object['chain_ids'] = ag.chain_ids
        mol_object['atom_type_unique'] = ag.atom_type_unique
        mol_object.mn['subframes'] = subframes
        mol_object.mn['molecule_type'] = 'md'

        # add the atomgroup to the session
        # the name of the atomgroup may be different from
        # the input name because of the uniqueness requirement
        # of the blender object name.
        # instead, the name generated by blender is used.
        #if mol_object.name != name:
        #    warnings.warn(
        #        "The name of the object is changed to {} because {} is already used.".format(
        #            mol_object.name, name
        #        )
        #    )

        #self.atom_reps[mol_object.name] = ag
        #self.universe_reps[mol_object.name] = {
        #    "universe": ag.universe,
        #    "frame_mapping": frame_mapping,
        #}
        #self.rep_names.append(mol_object.name)

        # for old import, the node tree is added all at once
        # in the end of in_memory
        if node_setup:
            nodes.create_starting_node_tree(
                object=mol_object,
                style=style,
            )

        if return_object:
            return mol_object



@persistent
def _rejuvenate_universe(scene):
    """
    Rejuvenate the session when the old Blend file is loaded.
    It will search through all the objects in the scene and
    find the ones that are molecules.
    It requires the pkl file to be in the default location and
    still exist.

    Warning:
    -------
    When a Blend file saved from another computer is loaded,
    the session will likely be lost.
    The Blend file also need to be opened from the same place
    (working directory) as when it is saved.
    """
    mol_objects = {}
    for object in bpy.data.objects:
        try:
            obj_type = object["type"]
            if obj_type == "molecule":
                mol_objects[object.name] = object
        except KeyError:
            pass

    if len(mol_objects) > 0:
        bpy.types.Scene.mda_session = MDAnalysisSession._rejuvenate(
            mol_objects)


@persistent
def _sync_universe(scene):
    """
    Sync the universe when the Blend file is saved.
    It will be saved as a .mda_session file in the
    same place as the Blend file).
    """
    if hasattr(bpy.types.Scene, "mda_session"):
        blender_save_loc = bpy.data.filepath
        if bpy.types.Scene.mda_session is not None:
            bpy.types.Scene.mda_session._dump(blender_save_loc)

    
