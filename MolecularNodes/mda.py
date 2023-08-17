import bpy
from bpy.app.handlers import persistent

import MDAnalysis as mda
import numpy as np
import warnings
import pickle
import uuid
import os
from typing import Union, List, Tuple, Dict, Optional

from . import data
from . import coll
from . import obj
from . import nodes


class AtomGroupInBlender:
    def __init__(self,
                 ag,
                 representation='vdw',
                 include_bonds=True,
                 world_scale=0.01):
        """
        AtomGroup in Blender.

        Parameters:
        ----------
        ag : MDAnalysis.AtomGroup
            The atomgroup to add in the scene.
        representation : str, optional
            The representation of the atoms (default: "vdw").
        include_bonds : bool, optional
            Whether to include bond information if available (default: True).
        world_scale : float, optional
            The scaling factor for the world coordinates (default: 0.01).
        """

        self.ag = ag
        self.include_bonds = include_bonds
        self.world_scale = world_scale
        self.representation = representation

    @property
    def n_atoms(self):
        return self.ag.n_atoms
    
    @property
    def representation(self):
        return self._representation
    
    @representation.setter
    def representation(self, representation):
        if representation not in ['vdw']:
            raise ValueError("Representation can only be 'vdw' at the moment.")
        self._representation = representation
    
    
    @staticmethod
    def bool_selection(ag, selection):
        return np.isin(ag.ix, ag.select_atoms(selection).ix).astype(bool)
    
    @property
    def positions(self):
        return self.ag.positions * self.world_scale
 
    @property
    def bonds(self):
        if hasattr(self.ag, "bonds") and self.include_bonds:
            index_map = {index: i for i, index in enumerate(self.ag.indices)}

            new_bonds = []
            for bond in self.ag.bonds.indices:
                try:
                    new_index = [index_map[y] for y in bond]
                    new_bonds.append(new_index)
                except KeyError:
                    # fragment - one of the atoms in the bonds was
                    # deleted by the selection, so we shouldn't
                    # pass this as a bond.
                    pass

            bonds = np.array(new_bonds)
        else:
            bonds = []
        return bonds

    @property
    def elements(self):
        try:
            elements = self.ag.elements.tolist()
        except:
            try:
                elements = [
                    mda.topology.guessers.guess_atom_element(x) for x in self.ag.atoms.names
                ]
            except:
                elements = ['X'] * self.ag.n_atoms
        return elements

    @property
    def atomic_number(self):
        return np.array(
            [data.elements.get(element,
                              data.elements.get('X'))
                            .get('atomic_number') for element in self.elements]
        )

    @property
    def vdw_radii(self):
        # pm to Angstrom
        return np.array(
            [data.elements.get(element,
                              data.elements.get('X'))
                        .get('vdw_radii') for element in self.elements]) * 0.01 * self.world_scale

    @property
    def res_id(self):
        return self.ag.resnums

    @property
    def res_name(self):
        return np.array(list(map(lambda x: x[0:3], self.ag.resnames)))

    @property
    def res_num(self):
        return np.array(
            [data.residues.get(res_name,
                              data.residues.get('UNK'))
                            .get('res_name_num') for res_name in self.res_name]
        )

    @property
    def b_factor(self):
        if hasattr(self.ag, "tempfactors"):
            return self.ag.tempfactors
        else:
            return np.zeros(self.ag.n_atoms)

    @property
    def chain_id(self):
        if hasattr(self.ag, "chainIDs"):
            return self.ag.chainIDs
        else:
            return np.zeros(self.ag.n_atoms)
    
    @property
    def chain_id_unique(self):
        return np.unique(self.chain_id)

    @property
    def chain_id_num(self):
        chain_id_unique, chain_id_index = np.unique(self.chain_id, return_inverse=True)
        return chain_id_index

    @property
    def atom_type(self):
        return self.ag.types

    @property
    def atom_type_unique(self):
        return np.unique(self.atom_type)
    
    @property
    def atom_type_num(self):
        atom_type_unique, atom_type_index = np.unique(self.atom_type, return_inverse=True)
        return atom_type_index
    
    @property
    def is_nucleic(self):
        return self.bool_selection(self.ag, "nucleic")
    
    @property
    def is_peptide(self):
        return self.bool_selection(self.ag, "protein")
    
    @property
    def is_backbone(self):
        return self.bool_selection(self.ag, "backbone or nucleicbackbone")

    @property
    def is_alpha_carbon(self):
        return self.bool_selection(self.ag, "name CA")

    @property
    def is_solvent(self):
        return self.bool_selection(self.ag, "name OW or name HW1 or name HW2")
    
    @property
    def attributes(self):
        return {
            "atomic_number": {
                "value": self.atomic_number,
                "type": "INT",
                "domain": "POINT",
            },
            "vdw_radii": {
                "value": self.vdw_radii,
                "type": "FLOAT",
                "domain": "POINT",
            },
            "res_id": {
                "value": self.res_id,
                "type": "INT",
                "domain": "POINT",
            },
            "res_name": {
                "value": self.res_num,
                "type": "INT",
                "domain": "POINT",
            },
            "b_factor": {
                "value": self.b_factor,
                "type": "FLOAT",
                "domain": "POINT",
            },
            "chain_id": {
                "value": self.chain_id_num,
                "type": "INT",
                "domain": "POINT",
            },
            "atom_types": {
                "value": self.atom_type_num,
                "type": "INT",
                "domain": "POINT",
            },
            "is_backbone": {
                "value": self.is_backbone,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            "is_alpha_carbon": {
                "value": self.is_alpha_carbon,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            "is_solvent": {
                "value": self.is_solvent,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            "is_nucleic": {
                "value": self.is_nucleic,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            "is_peptide": {
                "value": self.is_peptide,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
        }
    
    

    
class MDAnalysisSession:
    # default location to store the session files
    session_tmp_dir = f"{os.path.expanduser('~')}/.blender_mda_session/"

    def __init__(self, world_scale=0.01):
        """
        Initialize a MDAnalysisSession.

        A unique uuid is generated for each session.
        During saving, the session is saved in the
        default location (`~/.blender_mda_session/` as a pickle file.

        Parameters:
        ----------
        world_scale : float, optional
            The scaling factor for the world coordinates (default: 0.01).
        
            
        Attributes:
        ----------
        world_scale : float
            The scaling factor for the world coordinates.   
        uuid : str
            The unique identifier for the session.
        universe_reps : dict
            A dictionary of the universes in the session.
        atom_reps : dict
            A dictionary of the atom representations in the session.
        rep_names : list
            A list of the names of the representations in the session.

        """

        # if the session already exists, load the existing session
        if hasattr(bpy.types.Scene, "mda_session"):
            warnings.warn("The existing mda session is loaded.")
            existing_session = bpy.types.Scene.mda_session
            self.__dict__ = existing_session.__dict__
            return

        self.world_scale = world_scale
        os.makedirs(self.session_tmp_dir, exist_ok=True)
        self.uuid = str(uuid.uuid4().hex)
        self.universe_reps = {}
        self.atom_reps = {}
        self.rep_names = []

        bpy.types.Scene.mda_session = self

        bpy.app.handlers.frame_change_post.append(
            self.update_trajectory_handler_wrapper()
        )

    def show(
        self,
        atoms : Union[mda.Universe, mda.AtomGroup],
        representation : str = "vdw",
        selection : str = "all",
        name : str = "atoms",
        include_bonds : bool = True,
        custom_selections : Dict[str, str] = {},
        frame_offset : int = 0,
    ):
        """
        Display an `MDAnalysis.Universe` or
       `MDAnalysis.core.groups.Atomgroup` in Blender.

        Parameters:
        ----------
        atoms : MDAnalysis.Universe or MDAnalysis.core.groups.Atomgroup
            The universe to load into blender.
        representation : str, optional
            The representation of the atoms
            (default: "vdw").
        selection : str, optional
            The selection string for atom filtering
            (default: "all").
            Uses MDAnalysis selection syntax.
        name : str, optional
            The name of the default atoms
            (default: "atoms").
        include_bonds : bool, optional
            Whether to include bond information if available
            (default: True).
        custom_selections : dict, optional
            A dictionary of custom selections for atom filtering with
            {'name' : 'selection string'}
            (default: {}).
            Uses MDAnalysis selection syntax.
        frame_offset : int, optional
            The frame offset for the trajectory.
            It means the frame number in Blender will be
            the absolute frame number minus the frame_offset
            (default: 0).
        """

        if isinstance(atoms, mda.Universe):
            atoms = atoms.select_atoms(selection)
            
        universe = atoms.universe

        mol_object = self._process_atomgroup(
                    ag=atoms,
                    name=name,
                    representation=representation,
                    include_bonds=include_bonds,
                    frame_offset=frame_offset,
                    return_object=True)
        
        # add the custom selections if they exist
        for name, sel in custom_selections.items():
            try:
                ag = universe.select_atoms(sel)
                if ag.n_atoms == 0:
                    raise ValueError("Selection is empty")
                self._process_atomgroup(
                    ag=ag, 
                    name=name,
                    representation=representation,
                    include_bonds=include_bonds,
                    frame_offset=frame_offset,
                    return_object=False
                    )
            except ValueError:
                warnings.warn("Unable to add custom selection: {}".format(sel.name))

        # check if the universe is already in the session

        bpy.context.view_layer.objects.active = mol_object

    def _process_atomgroup(self,
                      ag,
                      name='atoms',
                      representation='vdw',
                      include_bonds=True,
                      frame_offset=0,
                      return_object=False):
        """
        add the atomgroup in the Blender scene.

        Parameters:
        ----------
        ag : MDAnalysis.AtomGroup
            The atomgroup to add in the scene.
        name : str
            The name of the atomgroup. Default: 'atoms'
        representation : str
            The representation of the atoms. Default: 'vdw'
        include_bonds : bool
            Whether to include bond information if available. Default: True
        frame_offset : int
            The frame offset for the trajectory. Default: 0
        return_object : bool
            Whether to return the blender object or not. Default: False
        """
        ag_blender = AtomGroupInBlender(
                                        ag=ag,
                                        include_bonds=include_bonds,
                                        representation=representation,
                                        world_scale=self.world_scale
                                        )
        # create the initial model
        mol_object = obj.create_object(
            name=name,
            collection=coll.mn(),
            locations=ag_blender.positions,
            bonds=ag_blender.bonds,
        )

        mol_object["session"] = self.uuid

        ## add the attributes for the model
        for name, att in ag_blender.attributes.items():
            obj.add_attribute(
                mol_object, name, att["value"], att["type"], att["domain"]
            )


        mol_object['chain_id_unique'] = ag_blender.chain_id_unique
        mol_object['atom_type_unique'] = ag_blender.atom_type_unique

        self.atom_reps[mol_object.name] = ag_blender
        self.universe_reps[mol_object.name] = {
            "universe": ag.universe,
            "frame_offset": frame_offset,
        }
        self.rep_names.append(mol_object.name)

        nodes.create_starting_node_tree(
            obj=mol_object, starting_style=bpy.context.scene.mol_import_default_style
        )

        if return_object:
            return mol_object

    @persistent
    def update_trajectory(self, frame):
        for name in self.rep_names:
            universe = self.universe_reps[name]["universe"]
            frame_offset = self.universe_reps[name]["frame_offset"]
            if frame - frame_offset < 0:
                continue
            if universe.trajectory.n_frames <= frame - frame_offset:
                continue

            # only load the frame if it's not already loaded
            if not universe.trajectory.frame == frame - frame_offset:
                universe.trajectory[frame - frame_offset]

            ag_rep = self.atom_reps[name]
            mol_object = bpy.data.objects[name]
            locations = ag_rep.positions

            # if the class of AtomGroup is UpdatingAtomGroup
            # then update as a new mol_object
            if isinstance(ag_rep.ag, mda.core.groups.UpdatingAtomGroup):
                mol_object.data.clear_geometry()
                mol_object.data.from_pydata(
                                    ag_rep.positions,
                                    ag_rep.bonds,
                                    faces=[])
                for name, att in ag_rep.attributes.items():
                    obj.add_attribute(
                        mol_object, name, att["value"], att["type"], att["domain"]
                    )


                mol_object['chain_id_unique'] = ag_rep.chain_id_unique
                mol_object['atom_type_unique'] = ag_rep.atom_type_unique
                # add the attributes for the model

            else:
                # The only gotcha is that currently to write to vector
                # attributes such as position,
                # you have to supply a 1D array so we have to just reshape
                # it before setting. Not doing this can segfault blender.
                # They have plans to improve this API eventually but this
                # is what we currently have to work with.

                # for vert, loc in zip(mol_object.data.vertices, locations):
                #     vert.co = loc
                mol_object.data.attributes['position'].data.foreach_set(
                                        'vector',
                                        locations.reshape(-1))
                mol_object.data.update()

    @persistent
    def update_trajectory_handler_wrapper(self):
        def update_trajectory_handler(scene):
            self.remove_deleted_mol_objects()
            frame = scene.frame_current
            self.update_trajectory(frame)

        return update_trajectory_handler

    @persistent
    def remove_deleted_mol_objects(self):
        for name in self.rep_names:
            if name not in bpy.data.objects:
                self.rep_names.remove(name)
                del self.atom_reps[name]
                del self.universe_reps[name]

    def dump(self):
        with open(f"{self.session_tmp_dir}/{self.uuid}.pkl", "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def rejuvenate(cls, mol_objects):
        # get session name from mol_objects dictionary
        session_name = mol_objects[list(mol_objects.keys())[0]]['session']
        with open(f"{cls.session_tmp_dir}/{session_name}.pkl", "rb") as f:
            cls = pickle.load(f)
        bpy.app.handlers.frame_change_post.append(
            cls.update_trajectory_handler_wrapper()
        )
        return cls


@persistent
def rejuvenate_universe(scene):
    mol_objects = {}
    for object in bpy.data.objects:
        try:
            obj_type = object["type"]
            if obj_type == "molecule":
                mol_objects[object.name] = object
        except KeyError:
            pass

    if len(mol_objects) > 0:
        bpy.types.Scene.mda_session = MDAnalysisSession.rejuvenate(mol_objects)


@persistent
def sync_universe(scene):
    if bpy.types.Scene.mda_session is not None:
        bpy.types.Scene.mda_session.dump()