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
                 ag: mda.AtomGroup,
                 representation: str = "vdw",
                 include_bonds: bool = True,
                 world_scale: float = 0.01):
        """
        AtomGroup in Blender.
        It will be dynamically updated when the frame changes or
        when the topology of the underlying atoms changes.

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

        Attributes:
        ----------
        ag : MDAnalysis.AtomGroup
            The atomgroup to add in the scene.
        include_bonds : bool
            Whether to include bond information if available.
        world_scale : float
            The scaling factor for the world coordinates.
        representation : str
            The representation of the atoms in Blender.
            Currently only 'vdw' is supported.
        n_atoms : int
            The number of atoms in the atomgroup.
        positions : np.ndarray
            The positions of the atoms in the atomgroup.
        bonds : list
            The bonds of the atoms in the atomgroup.
            If include_bonds is False, then bonds is an empty list.
        elements : list
            The elements of the atoms in the atomgroup.
            If the elements are not available,
            then the elements are guessed from the atom names.
        atomic_number : np.ndarray
            The atomic numbers of the atoms in the atomgroup
            based on their elements.
        vdw_radii : np.ndarray
            The van der Waals radii of the atoms in the atomgroup
            based on their elements.
        res_id : np.ndarray
            The residue ids of the atoms in the atomgroup.
        res_name : np.ndarray
            The residue names of the atoms in the atomgroup.
        res_num : np.ndarray
            The residue numbers of the atoms in the atomgroup
            The residue numbers are based on the residue names
            and are stored in the data.residues dictionary.
        b_factor : np.ndarray
            The B-factors of the atoms in the atomgroup.
        chain_id : np.ndarray
            The chain ids of the atoms in the atomgroup.
        chain_id_unique : np.ndarray
            The unique chain ids of the atoms in the atomgroup.
        chain_id_num : np.ndarray
            The chain id numbers of the atoms in the atomgroup.
            It is the index of the unique chain ids.
        atom_type : np.ndarray
            The atom types of the atoms in the atomgroup.
        atom_type_unique : np.ndarray
            The unique atom types of the atoms in the atomgroup.
        atom_type_num : np.ndarray
            The atom type numbers of the atoms in the atomgroup.
            It is the index of the unique atom types.
        is_nucleic : np.ndarray
            Whether the atoms in the atomgroup are nucleic.
        is_peptide : np.ndarray
            Whether the atoms in the atomgroup are peptide.
        is_backbone : np.ndarray
            Whether the atoms in the atomgroup are backbone.
        is_alpha_carbon : np.ndarray
            Whether the atoms in the atomgroup are alpha carbon.
        is_solvent : np.ndarray
            Whether the atoms in the atomgroup are solvent.
        """

        self.ag = ag
        self.include_bonds = include_bonds
        self.world_scale = world_scale
        self.representation = representation

    @property
    def n_atoms(self) -> int:
        return self.ag.n_atoms
    
    @property
    def representation(self) -> str:
        return self._representation
    
    @representation.setter
    def representation(self, representation):
        if representation not in ['vdw']:
            raise ValueError("Representation can only be 'vdw' at the moment.")
        self._representation = representation
    
    @staticmethod
    def bool_selection(ag, selection) -> np.ndarray:
        return np.isin(ag.ix, ag.select_atoms(selection).ix).astype(bool)
    
    @property
    def positions(self) -> np.ndarray:
        return self.ag.positions * self.world_scale
 
    @property
    def bonds(self) -> List[List[int]]:
        if hasattr(self.ag, "bonds") and self.include_bonds:
            bond_indices = self.ag.bonds.indices
            atm_indices = self.ag.indices
            bond_filtering = np.all(np.isin(bond_indices, atm_indices), axis=1)
            bond_indices = bond_indices[bond_filtering]

            index_map = {index: i for i, index in enumerate(self.ag.indices)}

            bonds = [[index_map[bond[0]], index_map[bond[1]]] for bond in bond_indices]
        else:
            bonds = []
        return bonds

    @property
    def elements(self) -> List[str]:
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
    def atomic_number(self) -> np.ndarray:
        return np.array(
            [data.elements.get(element,
                              data.elements.get('X'))
                            .get('atomic_number') for element in self.elements]
        )

    @property
    def vdw_radii(self) -> np.ndarray:
        # pm to Angstrom
        return np.array(
            [data.elements.get(element,
                              data.elements.get('X'))
                        .get('vdw_radii') for element in self.elements]) * 0.01 * self.world_scale

    @property
    def res_id(self) -> np.ndarray:
        return self.ag.resnums

    @property
    def res_name(self) -> np.ndarray:
        return np.array(list(map(lambda x: x[0:3], self.ag.resnames)))

    @property
    def res_num(self) -> np.ndarray:
        return np.array(
            [data.residues.get(res_name,
                              data.residues.get('UNK'))
                            .get('res_name_num') for res_name in self.res_name]
        )

    @property
    def b_factor(self) -> np.ndarray:
        if hasattr(self.ag, "tempfactors"):
            return self.ag.tempfactors
        else:
            return np.zeros(self.ag.n_atoms)

    @property
    def chain_id(self) -> np.ndarray:
        if hasattr(self.ag, "chainIDs"):
            return self.ag.chainIDs
        else:
            return np.zeros(self.ag.n_atoms)
    
    @property
    def chain_id_unique(self) -> np.ndarray:
        return np.unique(self.chain_id)

    @property
    def chain_id_num(self) -> np.ndarray:
        chain_id_unique, chain_id_index = np.unique(self.chain_id, return_inverse=True)
        return chain_id_index

    @property
    def atom_type(self) -> np.ndarray:
        return self.ag.types

    @property
    def atom_type_unique(self) -> np.ndarray:
        return np.unique(self.atom_type)
    
    @property
    def atom_type_num(self) -> np.ndarray:
        atom_type_unique, atom_type_index = np.unique(self.atom_type, return_inverse=True)
        return atom_type_index
    
    @property
    def is_nucleic(self) -> np.ndarray:
        return self.bool_selection(self.ag, "nucleic")
    
    @property
    def is_peptide(self) -> np.ndarray:
        return self.bool_selection(self.ag, "protein")
    
    @property
    def is_backbone(self) -> np.ndarray:
        return self.bool_selection(self.ag, "backbone or nucleicbackbone")

    @property
    def is_alpha_carbon(self) -> np.ndarray:
        return self.bool_selection(self.ag, "name CA")

    @property
    def is_solvent(self) -> np.ndarray:
        return self.bool_selection(self.ag, "name OW or name HW1 or name HW2")
    
    @property
    def _attributes_2_blender(self):
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
    """
    The MDAnalysis session.
    
    The MDAnalysis session is the main class that stores the
    MDAnalysis data in Blender.
    It is a singleton class that is initialized when the first
    MDAnalysis data is loaded in Blender.
    The MDAnalysis session is loaded when Blender is restarted.
    The MDAnalysis session is updated when the frame changes.
    The MDAnalysis session is dumped when a Blender file is saved.
    The MDAnalysis session is saved as a pickle file in the
    default location (`~/.blender_mda_session/`).

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
    session_tmp_dir : str
        The default location to store the session files.
    """

    # default location to store the session files
    session_tmp_dir = f"{os.path.expanduser('~')}/.blender_mda_session/"

    def __init__(self, world_scale: float = 0.01, legacy: bool = False):
        """
        Initialize a MDAnalysisSession.

        A unique uuid is generated for each session.
        During saving, the session is saved in the
        default location (`~/.blender_mda_session/` as a pickle file.

        Parameters:
        ----------
        world_scale : float, optional
            The scaling factor for the world coordinates (default: 0.01).
        legacy : bool, optional
            Whether the old import is used (default: False).
        """

        # if the session already exists, load the existing session
        if hasattr(bpy.types.Scene, "mda_session"):
            warnings.warn("The existing mda session is loaded.")
            existing_session = bpy.types.Scene.mda_session
            self.__dict__ = existing_session.__dict__
            return

        self.world_scale = world_scale
        self.universe_reps = {}
        self.atom_reps = {}
        self.rep_names = []
        self.uuid = str(uuid.uuid4().hex)

        if legacy:
            return
        os.makedirs(self.session_tmp_dir, exist_ok=True)
        bpy.types.Scene.mda_session = self
        bpy.app.handlers.frame_change_post.append(
            self._update_trajectory_handler_wrapper()
        )
        bpy.app.handlers.depsgraph_update_pre.append(
            self._update_representations_handler_wrapper()
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
       `MDAnalysis.Atomgroup` in Blender.

        Parameters:
        ----------
        atoms : MDAnalysis.Universe or MDAnalysis.Atomgroup
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

        if representation not in ['vdw']:
            warnings.warn("Representation can only be 'vdw' at the moment.")

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
        for sel_name, sel in custom_selections.items():
            try:
                ag = universe.select_atoms(sel)
                if ag.n_atoms == 0:
                    raise ValueError("Selection is empty")
                self._process_atomgroup(
                    ag=ag, 
                    name=sel_name,
                    representation=representation,
                    include_bonds=include_bonds,
                    frame_offset=frame_offset,
                    return_object=False
                    )
            except ValueError:
                warnings.warn("Unable to add custom selection: {}".format(name))

        bpy.context.view_layer.objects.active = mol_object

    def _process_atomgroup(
        self,
        ag,
        name="atoms",
        representation="vdw",
        include_bonds=True,
        frame_offset=0,
        add_node_tree=True,
        return_object=False,
    ):
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
        add_node_tree : bool
            Whether to add the node tree for the atomgroup. Default: True
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

        # add the attributes for the model in blender
        for att_name, att in ag_blender._attributes_2_blender.items():
            obj.add_attribute(
                mol_object, att_name, att["value"], att["type"], att["domain"]
            )
        mol_object['chain_id_unique'] = ag_blender.chain_id_unique
        mol_object['atom_type_unique'] = ag_blender.atom_type_unique

        # add the atomgroup to the session
        # the name of the atomgroup may be different from
        # the input name because of the uniqueness requirement
        # of the blender object name.
        # instead, the name generated by blender is used.
        if mol_object.name != name:
            warnings.warn(
                "The name of the object is changed to {} because {} is already used.".format(
                    mol_object.name, name
                )
            )

        self.atom_reps[mol_object.name] = ag_blender
        self.universe_reps[mol_object.name] = {
            "universe": ag.universe,
            "frame_offset": frame_offset,
        }
        self.rep_names.append(mol_object.name)

        if add_node_tree:
            nodes.create_starting_node_tree(
                obj=mol_object,
                starting_style=bpy.context.scene.mol_import_default_style,
            )

        if return_object:
            return mol_object

    @persistent
    def _update_trajectory(self, frame):
        for rep_name in self.rep_names:
            universe = self.universe_reps[rep_name]["universe"]
            frame_offset = self.universe_reps[rep_name]["frame_offset"]
            if frame - frame_offset < 0:
                continue
            if universe.trajectory.n_frames <= frame - frame_offset:
                continue

            # only load the frame if it's not already loaded
            if not universe.trajectory.frame == frame - frame_offset:
                universe.trajectory[frame - frame_offset]

            ag_rep = self.atom_reps[rep_name]
            mol_object = bpy.data.objects[rep_name]
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
                    obj.add_attribute(
                        mol_object, att_name, att["value"], att["type"], att["domain"]
                    )
                mol_object['chain_id_unique'] = ag_rep.chain_id_unique
                mol_object['atom_type_unique'] = ag_rep.atom_type_unique
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
    def _update_trajectory_handler_wrapper(self):
        def update_trajectory_handler(scene):
            frame = scene.frame_current
            self._update_trajectory(frame)

        return update_trajectory_handler

    @persistent
    def _update_representations_handler_wrapper(self):
        def update_representations_handler(scene):
            self._remove_deleted_mol_objects()
            #TODO: check for topology changes
            #TODO: update for representation changes
        
        return update_representations_handler

    @persistent
    def _remove_deleted_mol_objects(self):
        for rep_name in self.rep_names:
            if rep_name not in bpy.data.objects:
                self.rep_names.remove(rep_name)
                del self.atom_reps[rep_name]
                del self.universe_reps[rep_name]

    def show_legacy(
        self,
        atoms: Union[mda.Universe, mda.AtomGroup],
        representation: str = "vdw",
        selection: str = "all",
        name: str = "atoms",
        include_bonds: bool = True,
        custom_selections: Dict[str, str] = {},
    ):
        print(custom_selections)
        if isinstance(atoms, mda.Universe):
            atoms = atoms.select_atoms(selection)

        universe = atoms.universe

        mol_object = self._process_atomgroup(
            ag=atoms,
            name=name,
            representation=representation,
            include_bonds=include_bonds,
            add_node_tree=False,
            return_object=True,
        )

        for sel_name, sel in custom_selections.items():
            obj.add_attribute(
                object=mol_object,
                name=sel_name,
                data=AtomGroupInBlender.bool_selection(atoms, sel),
                type="BOOLEAN",
                domain="POINT",
            )

        coll_frames = coll.frames(name)

        add_occupancy = True
        for ts in universe.trajectory:
            frame = obj.create_object(
                name=name + "_frame_" + str(ts.frame),
                collection=coll_frames,
                locations=atoms.positions * self.world_scale,
            )
            # adds occupancy data to each frame if it exists
            # This is mostly for people who want to store frame-specific information in the
            # b_factor but currently neither biotite nor MDAnalysis give access to frame-specific
            # b_factor information. MDAnalysis gives frame-specific access to the `occupancy`
            # so currently this is the only method to get frame-specific data into MN
            # for more details: https://github.com/BradyAJohnston/MolecularNodes/issues/128
            if add_occupancy:
                try:
                    obj.add_attribute(frame, "occupancy", ts.data["occupancy"])
                except:
                    add_occupancy = False
        # disable the frames collection from the viewer
        bpy.context.view_layer.layer_collection.children[coll.mn().name].children[
            coll_frames.name
        ].exclude = True

        nodes.create_starting_node_tree(
            obj=mol_object,
            coll_frames=coll_frames,
            starting_style=bpy.context.scene.mol_import_default_style,
        )

        bpy.context.view_layer.objects.active = mol_object

    def transfer_to_memory(
        self, start=None, stop=None, step=None, verbose=False, **kwargs
    ):
        """
        Transfer the trajectories in the session to memory.

        Parameters:
        ----------
        start : int, optional
            The first frame to transfer (default: None).
            If None, then the first frame of the trajectory is used.
        stop : int, optional
            The last frame to transfer (default: None).
            If None, then the last frame of the trajectory is used.
        step : int, optional
            The step between frames (default: None).
            If None, then the step is 1.
        verbose : bool, optional
            Whether to print the progress (default: False).
        """

        warnings.warn("The trajectories in this session \n"
                      "is transferred to memory. \n"
                      "All the frame information will be saved in \n"
                      "the tmp file ~/.blender_mda_session/ \n"
                      f"{self.uuid}.pkl when the blend file \n"
                      "is saved.")

        for rep_name in self.rep_names:
            universe = self.universe_reps[rep_name]["universe"]
            universe.transfer_to_memory(
                start=start, stop=stop, step=step, verbose=verbose, **kwargs
            )

    def _dump(self):
        with open(f"{self.session_tmp_dir}/{self.uuid}.pkl", "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def _rejuvenate(cls, mol_objects):
        # get session name from mol_objects dictionary
        session_name = mol_objects[list(mol_objects.keys())[0]]['session']
        with open(f"{cls.session_tmp_dir}/{session_name}.pkl", "rb") as f:
            cls = pickle.load(f)
        bpy.app.handlers.frame_change_post.append(
            cls._update_trajectory_handler_wrapper()
        )
        bpy.app.handlers.depsgraph_update_pre.append(
            cls._update_representations_handler_wrapper()
        )
        return cls


@persistent
def _rejuvenate_universe(scene):
    mol_objects = {}
    for object in bpy.data.objects:
        try:
            obj_type = object["type"]
            if obj_type == "molecule":
                mol_objects[object.name] = object
        except KeyError:
            pass

    if len(mol_objects) > 0:
        bpy.types.Scene.mda_session = MDAnalysisSession._rejuvenate(mol_objects)


@persistent
def _sync_universe(scene):
    if bpy.types.Scene.mda_session is not None:
        bpy.types.Scene.mda_session._dump()
