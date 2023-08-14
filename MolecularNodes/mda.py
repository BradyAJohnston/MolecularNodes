import bpy
from bpy.app.handlers import persistent

import MDAnalysis as mda
import numpy as np
import warnings
import pickle
import uuid
import os

from . import data
from . import coll
from . import obj
from . import nodes


class MDAnalysisSession:
    # default location to store the session files
    session_tmp_dir = f"{os.path.expanduser('~')}/.blender_mda_session/"

    def __init__(
        self,
        universe,
        name="atoms",
        md_start=0,
        md_end=-1,
        md_step=1,
        world_scale=0.01,
        include_bonds=True,
        selection="not (name H* or name OW)",
        custom_selections=None,
    ):
        """
        MDAnalysis.universe trajectory in Blender.

        Parameters:
        ----------
        universe : MDAnalysis.Universe
            The universe to load into blender.
        name : str, optional
            The name of the default atoms (default: "atoms").
        md_start : int, optional
            The starting frame of the trajectory to load (default: 0).
        md_end : int, optional
            The ending frame of the trajectory to load (default: 49).
        md_step : int, optional
            The step size between frames to load (default: 1).
        world_scale : float, optional
            The scaling factor for the world coordinates (default: 0.01).
        include_bonds : bool, optional
            Whether to include bond information if available (default: True).
        selection : str, optional
            The selection string for atom filtering (default: "not (name H* or name OW)").
            Uses MDAnalysis selection syntax.
        custom_selections : dict or None, optional
            A dictionary of custom selections for atom filtering with
            {'name' : 'selection string'} (default: None).

        Raises:
        ------
        FileNotFoundError
            If the topology or trajectory file is not found.
        IOError
            If there is an error reading the files.
        """
        self.universe = universe
        self.name = name
        self.start = md_start
        self.end = md_end
        self.step = md_step
        self.world_scale = world_scale
        self.include_bonds = include_bonds
        self.selection = selection
        os.makedirs(self.session_tmp_dir, exist_ok=True)
        #        self.custom_selections = custom_selections

        self.uuid = str(uuid.uuid4().hex)

        self.ags = {}
        self.mol_objects = {}
        self.ag_names = []
        self.load_universe()

        ag = self.universe.select_atoms(self.selection)
        self.add_atomgroup(ag=ag,
                           name=name)

        # add the custom selections if they exist
        for name, sel in custom_selections.items():
            try:
                ag = self.universe.select_atoms(sel)
                self.add_atomgroup(
                    ag=ag, 
                    name=name
                    )
            except:
                warnings.warn("Unable to add custom selection: {}".format(sel.name))

        bpy.app.handlers.frame_change_post.append(
            self.update_trajectory_handler_wrapper()
        )
        bpy.types.Scene.mda_session = self
        # store the universe, atomgroup, trajectory in the scene
        # so that they can be accessed later in Blender
        # TODO what if there are multiple universes?
        #        bpy.types.Scene.universe = self


    def load_universe(self):
        # separate the trajectory, separate to the topology or the subsequence selections
        self.trajectory = self.universe.trajectory[self.start : self.end : self.step]
        self.n_frames = len(self.trajectory)

    def add_atomgroup(self, ag, name='atoms'):
        """
        add the atomgroup in the Blender scene.

        Parameters:
        ----------
        ag : MDAnalysis.AtomGroup
            The atomgroup to add in the scene.
        name : str
            The name of the atomgroup. Default: 'atoms'
        """

        # Try and extract the elements from the topology. If the universe doesn't contain
        # the element information, then guess based on the atom names in the topology
        try:
            elements = ag.elements.tolist()
        except:
            try:
                elements = [
                    mda.topology.guessers.guess_atom_element(x) for x in ag.atoms.names
                ]
            except:
                pass

        if hasattr(ag, "bonds") and self.include_bonds:
            # If there is a selection, we need to recalculate the bond indices
            if self.selection != "":
                index_map = {index: i for i, index in enumerate(ag.indices)}

                new_bonds = []
                for bond in ag.bonds.indices:
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
                bonds = ag.bonds.indices

        else:
            bonds = []

        # create the initial model
        mol_object = obj.create_object(
            name=name,
            collection=coll.mn(),
            locations=ag.positions * self.world_scale,
            bonds=bonds,
        )

        mol_object["session"] = self.uuid

        ## add the attributes for the model

        # The attributes for the model are initially defined as single-use functions. This allows
        # for a loop that attempts to add each attibute by calling the function. Only during this
        # loop will the call fail if the attribute isn't accessible, and the warning is reported
        # there rather than setting up a try: except: for each individual attribute which makes
        # some really messy code.

        def att_atomic_number():
            atomic_number = np.array(
                list(
                    map(
                        # if getting the element fails for some reason, return an atomic number of -1
                        lambda x: data.elements.get(x, {"atomic_number": -1}).get(
                            "atomic_number"
                        ),
                        np.char.title(elements),
                    )
                )
            )
            return atomic_number

        def att_vdw_radii():
            try:
                vdw_radii = np.array(
                    list(
                        map(
                            lambda x: mda.topology.tables.vdwradii.get(x, 1),
                            np.char.upper(elements),
                        )
                    )
                )
            except:
                # if fail to get radii, just return radii of 1 for everything as a backup
                vdw_radii = np.ones(len(ag.names))
                warnings.warn(
                    "Unable to extract VDW Radii. Defaulting to 1 for all points."
                )

            return vdw_radii * self.world_scale

        def att_res_id():
            return ag.resnums

        def att_res_name():
            res_names = np.array(list(map(lambda x: x[0:3], ag.resnames)))
            res_numbers = np.array(
                list(
                    map(
                        lambda x: data.residues.get(x, {"res_name_num": 0}).get(
                            "res_name_num"
                        ),
                        res_names,
                    )
                )
            )
            return res_numbers

        def att_b_factor():
            return ag.tempfactors

        def att_chain_id():
            chain_id = ag.chainIDs
            chain_id_unique = np.unique(chain_id)
            chain_id_num = np.array(
                list(map(lambda x: np.where(x == chain_id_unique)[0][0], chain_id))
            )
            mol_object["chain_id_unique"] = chain_id_unique
            return chain_id_num

        # returns a numpy array of booleans for each atom, whether or not they are in that selection
        def bool_selection(selection):
            return np.isin(ag.ix, ag.select_atoms(selection).ix).astype(bool)

        def att_is_backbone():
            return bool_selection("backbone or nucleicbackbone")

        def att_is_alpha_carbon():
            return bool_selection("name CA")

        def att_is_solvent():
            return bool_selection("name OW or name HW1 or name HW2")

        def att_atom_type():
            return np.array(ag.types, dtype=int)

        def att_is_nucleic():
            return bool_selection("nucleic")

        def att_is_peptide():
            return bool_selection("protein")

        attributes = (
            {
                "name": "atomic_number",
                "value": att_atomic_number,
                "type": "INT",
                "domain": "POINT",
            },
            {
                "name": "vdw_radii",
                "value": att_vdw_radii,
                "type": "FLOAT",
                "domain": "POINT",
            },
            {"name": "res_id", "value": att_res_id, "type": "INT", "domain": "POINT"},
            {
                "name": "res_name",
                "value": att_res_name,
                "type": "INT",
                "domain": "POINT",
            },
            {
                "name": "b_factor",
                "value": att_b_factor,
                "type": "FLOAT",
                "domain": "POINT",
            },
            {
                "name": "chain_id",
                "value": att_chain_id,
                "type": "INT",
                "domain": "POINT",
            },
            {
                "name": "atom_types",
                "value": att_atom_type,
                "type": "INT",
                "domain": "POINT",
            },
            {
                "name": "is_backbone",
                "value": att_is_backbone,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            {
                "name": "is_alpha_carbon",
                "value": att_is_alpha_carbon,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            {
                "name": "is_solvent",
                "value": att_is_solvent,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            {
                "name": "is_nucleic",
                "value": att_is_nucleic,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            {
                "name": "is_peptide",
                "value": att_is_peptide,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
        )

        for att in attributes:
            # tries to add the attribute to the mesh by calling the 'value' function which returns
            # the required values do be added to the domain.
            try:
                obj.add_attribute(
                    mol_object, att["name"], att["value"](), att["type"], att["domain"]
                )
            except:
                warnings.warn(f"Unable to add attribute: {att['name']}.")

        self.ags[mol_object.name] = ag
        self.mol_objects[mol_object.name] = mol_object
        self.ag_names.append(mol_object.name)

        nodes.create_starting_node_tree(
            obj=mol_object, starting_style=bpy.context.scene.mol_import_default_style
        )

    @persistent
    def update_trajectory(self, frame):
        # TODO: how about a dynamic selection with variable number of atoms?
        self.trajectory[frame]
        for name in self.ag_names:
            ag = self.ags[name]
            mol_object = self.mol_objects[name]
            locations = ag.positions * self.world_scale
            for vert, loc in zip(mol_object.data.vertices, locations):
                vert.co = loc

    @persistent
    def update_trajectory_handler_wrapper(self):
        def update_trajectory_handler(scene):
            frame = scene.frame_current
            if frame >= self.n_frames:
                return None
            self.update_trajectory(frame)

        return update_trajectory_handler

    def __getstate__(self):
        state = self.__dict__.copy()
        # remove the unpicklable entries.
        del state["mol_objects"]
        return state
    
    def dump(self):
        with open(f"{self.session_tmp_dir}/{self.uuid}.pkl", "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def rejuvenate(cls, mol_objects):
        # get any object from mol_objects dictionary
        session_name = mol_objects[list(mol_objects.keys())[0]]['session']
        with open(f"{cls.session_tmp_dir}/{session_name}.pkl", "rb") as f:
            cls = pickle.load(f)
        cls.mol_objects = mol_objects
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

    # TODO for now we assume there is only one Universe
    if len(mol_objects) > 0:
        bpy.types.Scene.mda_session = MDAnalysisSession.rejuvenate(mol_objects)


@persistent
def sync_universe(scene):
    if bpy.types.Scene.mda_session is not None:
        bpy.types.Scene.mda_session.dump()