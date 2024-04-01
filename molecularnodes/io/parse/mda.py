import bpy
from bpy.app.handlers import persistent

try:
    import MDAnalysis as mda
except ImportError:
    HAS_mda = False
    import types

    class MockAtomGroup:
        pass

    class MockUniverse:
        pass

    mda = types.ModuleType("MDAnalysis")
    mda.Universe = MockUniverse
    mda.AtomGroup = MockAtomGroup
    mda.core = types.ModuleType("core")
    mda.topology = types.ModuleType("topology")

else:
    HAS_mda = True
import numpy as np
import warnings
import pickle
from typing import Union, List, Dict

from ... import data
from ...pkg import start_logging
from ...blender import (
    coll, obj, nodes
)
from ...utils import lerp


class AtomGroupInBlender:
    def __init__(self,
                 ag: mda.AtomGroup,
                 style: str = "vdw",
                 world_scale: float = 0.01):
        """
        AtomGroup in Blender.
        It will be dynamically updated when the frame changes or
        when the topology of the underlying atoms changes.

        Parameters:
        ----------
        ag : MDAnalysis.AtomGroup
            The atomgroup to add in the scene.
        style : str, optional
            The style of the atoms (default: "vdw").
        world_scale : float, optional
            The scaling factor for the world coordinates (default: 0.01).

        Attributes:
        ----------
        ag : MDAnalysis.AtomGroup
            The atomgroup to add in the scene.
        world_scale : float
            The scaling factor for the world coordinates.
        style : str
            The style of the atoms in Blender.
        n_atoms : int
            The number of atoms in the atomgroup.
        positions : np.ndarray
            The positions of the atoms in the atomgroup.
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
        if not HAS_mda:
            raise ImportError("MDAnalysis is not installed.")
        self.ag = ag
        self.world_scale = world_scale
        self.style = style

    @property
    def n_atoms(self) -> int:
        return self.ag.n_atoms

    @property
    def style(self) -> str:
        return self._style

    @style.setter
    def style(self, style):
        self._style = style

    @staticmethod
    def bool_selection(ag, selection) -> np.ndarray:
        return np.isin(ag.ix, ag.select_atoms(selection).ix).astype(bool)

    @property
    def positions(self) -> np.ndarray:
        return self.ag.positions * self.world_scale

    @property
    def bonds(self) -> List[List[int]]:
        if hasattr(self.ag, "bonds"):
            bond_indices = self.ag.bonds.indices
            atm_indices = self.ag.indices
            bond_filtering = np.all(np.isin(bond_indices, atm_indices), axis=1)
            bond_indices = bond_indices[bond_filtering]

            index_map = {index: i for i, index in enumerate(self.ag.indices)}

            bonds = [[index_map[bond[0]], index_map[bond[1]]]
                     for bond in bond_indices]
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
                        x if x in data.elements.keys() else
                        mda.topology.guessers.guess_atom_element(x) for x in self.ag.atoms.names]

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
                               {'vdw_radii': 100})
             .get('vdw_radii') for element in self.elements]) * 0.01 * self.world_scale

    @property
    def mass(self) -> np.ndarray:
        # units: daltons
        try: 
            masses = np.array([x.mass for x in self.ag.atoms])
        except mda.exceptions.NoDataError:
            masses = np.array(
                    [data.elements.get(element,
                                       {'standard_mass': 0})
                     .get('standard_mass') for element in self.elements])
        return masses 

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
    def chain_ids(self) -> np.ndarray:
        return np.unique(self.chain_id)

    @property
    def chain_id_num(self) -> np.ndarray:
        chain_ids, chain_id_index = np.unique(
            self.chain_id, return_inverse=True)
        return chain_id_index

    @property
    def atom_type(self) -> np.ndarray:
        return self.ag.types

    @property
    def atom_type_unique(self) -> np.ndarray:
        return np.unique(self.atom_type)

    @property
    def atom_type_num(self) -> np.ndarray:
        atom_type_unique, atom_type_index = np.unique(
            self.atom_type, return_inverse=True)
        return atom_type_index

    @property
    def atom_name(self) -> np.ndarray:
        if hasattr(self.ag, "names"):
            return self.ag.names
        else:
            return np.zeros(self.ag.n_atoms)

    @property
    def atom_name_num(self) -> np.ndarray:
        if hasattr(self.ag, "names"):
            return np.array(list(map(lambda x: data.atom_names.get(x, -1), self.atom_name)))
        else:
            return np.repeat(-1, self.ag.n_atoms)

    @property
    def is_nucleic(self) -> np.ndarray:
        return self.bool_selection(self.ag, "nucleic")

    @property
    def is_peptide(self) -> np.ndarray:
        return self.bool_selection(self.ag, "protein or (name BB SC*)")

    @property
    def is_lipid(self) -> np.ndarray:
        return np.isin(self.ag.resnames, data.lipid_names)

    @property
    def is_backbone(self) -> np.ndarray:
        return self.bool_selection(self.ag, "backbone or nucleicbackbone or name BB")

    @property
    def is_alpha_carbon(self) -> np.ndarray:
        return self.bool_selection(self.ag, "name CA or name BB")

    @property
    def is_solvent(self) -> np.ndarray:
        return self.bool_selection(self.ag, "name OW or name HW1 or name HW2 or resname W or resname PW")

    @property
    def _attributes_2_blender(self):
        """
        The attributes that will be added to the Blender object.
        """
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
            "mass": {
                "value": self.mass,
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
            "atom_name": {
                "value": self.atom_name_num,
                "type": "INT",
                "domain": "POINT"
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
            "is_lipid": {
                "value": self.is_lipid,
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
    When a Blender file is saved, the MDAnalysis session will be
    dumped (with pickle) to the same place as the Blender file
    with .mda_session extension.

    Parameters:
    ----------
    world_scale : float, optional
        The scaling factor for the world coordinates (default: 0.01).

    Attributes:
    ----------
    world_scale : float
        The scaling factor for the world coordinates.
    universe_reps : dict
        A dictionary of the universes in the session.
    atom_reps : dict
        A dictionary of the atom styles in the session.
    rep_names : list
        A list of the names of the styles in the session.

    Methods:
    -------
    show(atoms, style, selection, name, custom_selections, frame_offset)
        Display an `MDAnalysis.Universe` or `MDAnalysis.Atomgroup` in Blender.
    in_memory(atoms, style, selection, name, custom_selections)
        Display an `MDAnalysis.Universe` or `MDAnalysis.Atomgroup` in Blender by loading all the
        frames as individual objects. Animation depends on the machinery inside geometric node.
    transfer_to_memory(start, stop, step, verbose, **kwargs)
        Transfer the trajectories in the session to memory.
    """

    def __init__(self, world_scale: float = 0.01, in_memory: bool = False):
        """
        Initialize a MDAnalysisSession.

        During saving, the session is pickled/serialized to the same
        location as the blend file with the extension .mda_session.
        The session is loaded when Blender is restarted.

        #TODO: Is it possible to start blender only when
        #a session is initialized? (Probably not for now)

        Parameters:
        ----------
        world_scale : float, optional
            The scaling factor for the world coordinates (default: 0.01).
        memory : bool, optional
            Whether the old import is used (default: False).
        """
        log = start_logging(logfile_name="mda")
        if not HAS_mda:
            raise ImportError("MDAnalysis is not installed.")

        # if the session already exists, load the existing session
        if hasattr(bpy.types.Scene, "mda_session"):
            warnings.warn("The existing mda session is loaded.")
            log.warning("The existing mda session is loaded.")
            existing_session = bpy.types.Scene.mda_session
            self.__dict__ = existing_session.__dict__
            return

        self.world_scale = world_scale
        self.universe_reps = {}
        self.atom_reps = {}
        self.rep_names = []

        if in_memory:
            return
        bpy.types.Scene.mda_session = self
        bpy.app.handlers.frame_change_post.append(
            self._update_trajectory_handler_wrapper()
        )
        bpy.app.handlers.depsgraph_update_pre.append(
            self._update_style_handler_wrapper()
        )
        log.info("MDAnalysis session is initialized.")

    @property
    def universe(self) -> mda.Universe:
        """
        The universe of the current active object.
        If the current active object is not an atom style,
        then the first atom style is used.
        """
        name = bpy.context.view_layer.objects.active.name
        try:
            return self.universe_reps[name]["universe"]
        except KeyError:
            return self.universe_reps[self.rep_names[0]]["universe"]

    def show(
        self,
        atoms: Union[mda.Universe, mda.AtomGroup],
        style: str = "vdw",
        selection: str = "all",
        name: str = "atoms",
        custom_selections: Dict[str, str] = {},
        frame_mapping: np.ndarray = None,
        subframes: int = 0,
        in_memory: bool = False
    ):
        """
        Display an `MDAnalysis.Universe` or
       `MDAnalysis.Atomgroup` in Blender.

        Parameters:
        ----------
        atoms : MDAnalysis.Universe or MDAnalysis.Atomgroup
            The universe to load into blender.
        style : str, optional
            The style to represent the atoms inside of Blender
            (default: "vdw").
        selection : str, optional
            The selection string for atom filtering
            (default: "all").
            Uses MDAnalysis selection syntax.
        name : str, optional
            The name of the default atoms
            (default: "atoms").
        custom_selections : dict, optional
            A dictionary of custom selections for atom filtering with
            {'name' : 'selection string'}
            (default: {}).
            Uses MDAnalysis selection syntax.
        frame_mapping : np.ndarray, optional
            A mapping from the frame indices in the Blender frame indices.
            for example a frame_mapping of [0, 0, 1, 1, 2, 3] will map
            the 1st frame (index 0) in the trajectory to the 1st and 2nd frames
            in Blender and so on.
            Note a subframes other than 1 will expand the frame_mapping from its
            original length to (subframes + 1) * original length.
            (default: None) which will map the frames in the trajectory
            to the frames in Blender one-to-one.
        subframes : int, optional
            The number of subframes to interpolate between each frame.
            (default: 0).
        in_memory : bool, optional
            Whether load the display in Blender by loading all the
            frames as individual objects.
            (default: False)
        """
        log = start_logging(logfile_name="mda")
        if in_memory:
            mol_object = self.in_memory(
                atoms=atoms,
                style=style,
                selection=selection,
                name=name,
                custom_selections=custom_selections
            )
            if frame_mapping is not None:
                warnings.warn("Custom frame_mapping not supported"
                              "when in_memory is on.")
            if subframes != 0:
                warnings.warn("Custom subframes not supported"
                              "when in_memory is on.")
            log.info(f"{atoms} is loaded in memory.")
            return mol_object
        if isinstance(atoms, mda.Universe):
            atoms = atoms.select_atoms(selection)

        universe = atoms.universe

        # if any frame_mapping is out of range, then raise an error
        if frame_mapping and (len(frame_mapping) > universe.trajectory.n_frames):
            raise ValueError("one or more mapping values are"
                             "out of range for the trajectory")

        mol_object = self._process_atomgroup(
            ag=atoms,
            frame_mapping=frame_mapping,
            subframes=subframes,
            name=name,
            style=style,
            return_object=True)

        # add the custom selections if they exist
        for sel_name, sel in custom_selections.items():
            try:
                ag = universe.select_atoms(sel)
                if ag.n_atoms == 0:
                    raise ValueError("Selection is empty")
                self._process_atomgroup(
                    ag=ag,
                    frame_mapping=frame_mapping,
                    subframes=subframes,
                    name=sel_name,
                    style=style,
                    return_object=False
                )
            except ValueError:
                warnings.warn(
                    "Unable to add custom selection: {}".format(name))

        bpy.context.view_layer.objects.active = mol_object
        log.info(f"{atoms} is loaded.")
        return mol_object

    def in_memory(
        self,
        atoms: Union[mda.Universe, mda.AtomGroup],
        style: str = "vdw",
        selection: str = "all",
        name: str = "atoms",
        custom_selections: Dict[str, str] = {},
        node_setup: bool = True
    ):
        """
        Display an `MDAnalysis.Universe` or
        `MDAnalysis.Atomgroup` in Blender by loading all the
        frames as individual objects. Animation depends on the machinery inside geometric node.

        Parameters:
        ----------
        atoms : MDAnalysis.Universe or MDAnalysis.Atomgroup
            The universe to load into blender.
        style : str, optional
            The style used to represent the atoms inside of Blender
            (default: "vdw").
        selection : str, optional
            The selection string for atom filtering
            (default: "all").
            Uses MDAnalysis selection syntax.
        name : str, optional
            The name of the default atoms
            (default: "atoms").
        custom_selections : dict, optional
            A dictionary of custom selections for atom filtering with
            {'name' : 'selection string'}
            (default: {}).
            Uses MDAnalysis selection syntax.
        node_setup : bool
            Whether to add the node tree for the atomgroup. Default: True
        """
        if isinstance(atoms, mda.Universe):
            atoms = atoms.select_atoms(selection)

        universe = atoms.universe

        mol_object = self._process_atomgroup(
            ag=atoms,
            name=name,
            style=style,
            node_setup=False,
            return_object=True,
        )

        for sel_name, sel in custom_selections.items():
            obj.set_attribute(
                object=mol_object,
                name=sel_name,
                data=AtomGroupInBlender.bool_selection(atoms, sel),
                type="BOOLEAN",
                domain="POINT",
            )

        coll_frames = coll.frames(name)

        # TODO: refractor it as a general feature
        add_occupancy = True
        for ts in universe.trajectory:
            frame = obj.create_object(
                name=name + "_frame_" + str(ts.frame),
                collection=coll_frames,
                vertices=atoms.positions * self.world_scale,
            )
            # adds occupancy data to each frame if it exists
            # This is mostly for people who want to store frame-specific information in the
            # b_factor but currently neither biotite nor MDAnalysis give access to frame-specific
            # b_factor information. MDAnalysis gives frame-specific access to the `occupancy`
            # so currently this is the only method to get frame-specific data into MN
            # for more details: https://github.com/BradyAJohnston/MolecularNodes/issues/128
            if add_occupancy:
                try:
                    obj.set_attribute(frame, "occupancy", ts.data["occupancy"])
                except:
                    add_occupancy = False

        # disable the frames collection from the viewer
        bpy.context.view_layer.layer_collection.children[coll.mn().name].children[
            coll_frames.name
        ].exclude = True

        if node_setup:
            nodes.create_starting_node_tree(
                object=mol_object,
                coll_frames=coll_frames,
                style=style,
            )

        bpy.context.view_layer.objects.active = mol_object

        return mol_object

    def transfer_to_memory(
        self, start=None, stop=None, step=None, verbose=False, **kwargs
    ):
        """
        Transfer the trajectories in the session to memory.
        This is an alternative way to make sure the blender session is
        independent of the original trajectory file.

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
        log = start_logging(logfile_name="mda")
        warnings.warn(
            "The trajectories in this session \n"
            "is transferred to memory. \n"
            "It is different from the in_memory loading \n"
            "because it uses the MemoryReader in MDAnalysis. \n"
            "instead of loading all the frames as individual objects. \n"
            "in Blender. \n"
        )

        for rep_name in self.rep_names:
            universe = self.universe_reps[rep_name]["universe"]
            universe.transfer_to_memory(
                start=start, stop=stop, step=step, verbose=verbose, **kwargs
            )
        log.info("The trajectories in this session is transferred to memory.")

    def _process_atomgroup(
        self,
        ag,
        frame_mapping=None,
        subframes=0,
        name="atoms",
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
        ag_blender = AtomGroupInBlender(
            ag=ag,
            style=style,
            world_scale=self.world_scale
        )
        # create the initial model
        mol_object = obj.create_object(
            name=name,
            collection=coll.mn(),
            vertices=ag_blender.positions,
            edges=ag_blender.bonds,
        )

        # add the attributes for the model in blender
        for att_name, att in ag_blender._attributes_2_blender.items():
            obj.set_attribute(
                mol_object, att_name, att["value"], att["type"], att["domain"]
            )
        mol_object['chain_ids'] = ag_blender.chain_ids
        mol_object['atom_type_unique'] = ag_blender.atom_type_unique
        mol_object.mn['subframes'] = subframes
        mol_object.mn['molecule_type'] = 'md'

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
            "frame_mapping": frame_mapping,
        }
        self.rep_names.append(mol_object.name)

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
