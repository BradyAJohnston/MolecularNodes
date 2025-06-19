import os
import pickle as pk
from typing import Dict, Union
import bpy
import MDAnalysis as mda
from bpy.app.handlers import persistent  # type: ignore
from bpy.props import StringProperty  # type: ignore
from bpy.types import Context  # type: ignore
from databpy.object import LinkedObjectError, get_from_uuid
from MDAnalysis.core.groups import AtomGroup
from .entities import Molecule
from .entities.ensemble.base import Ensemble
from .entities.trajectory.base import Trajectory
from .nodes.nodes import styles_mapping


def trim(dictionary: dict):
    dic = dictionary.copy()
    for key in list(dic.keys()):
        if dic[key].object is None:
            dic.pop(key)
    return dic


def make_paths_relative(trajectories: Dict[str, Trajectory]) -> None:
    for key, traj in trajectories.items():
        # save linked universe frame
        uframe = traj.uframe
        traj.universe.load_new(make_path_relative(traj.universe.trajectory.filename))
        # restore linked universe frame
        traj.uframe = uframe
        traj.save_filepaths_on_object()


def trim_root_folder(filename):
    "Remove one of the prefix folders from a filepath"
    return os.sep.join(filename.split(os.sep)[1:])


def make_path_relative(filepath):
    "Take a path and make it relative, in an actually usable way"
    try:
        filepath = os.path.relpath(filepath)
    except ValueError:
        return filepath

    # count the number of "../../../" there are to remove
    n_to_remove = int(filepath.count("..") - 2)
    # get the filepath without the huge number of "../../../../" at the start
    sans_relative = filepath.split("..")[-1]

    if n_to_remove < 1:
        return filepath

    for i in range(n_to_remove):
        sans_relative = trim_root_folder(sans_relative)

    return f"./{sans_relative}"


class MNSession:
    def __init__(self) -> None:
        self.entities: Dict[str, Union[Molecule, Trajectory, Ensemble]] = {}

    @property
    def molecules(self) -> Dict[str, Molecule]:
        return {k: v for k, v in self.entities.items() if isinstance(v, Molecule)}

    @property
    def trajectories(self) -> Dict[str, Trajectory]:
        # return a filtered dictionary of only the trajectories using isinstance(item, Trajectory)
        return {k: v for k, v in self.entities.items() if isinstance(v, Trajectory)}

    @property
    def ensembles(self) -> Dict[str, Ensemble]:
        # return a filtered dictionary of only the ensembles using isinstance(item, Ensemble)
        return {k: v for k, v in self.entities.items() if isinstance(v, Ensemble)}

    def register_entity(self, item: Union[Molecule, Trajectory, Ensemble]) -> None:
        """Add entity to the dictionary"""
        self.entities[item.uuid] = item

    def remove_entity(self, uuid: str) -> None:
        """Remove entity from the dictionary"""
        del self.entities[uuid]

    def match(self, obj: bpy.types.Object) -> Union[Molecule, Trajectory, Ensemble]:
        return self.get(obj.uuid)

    def get_object(self, uuid: str) -> bpy.types.Object | None:
        """
        Try and get an object from Blender's object database that matches the uuid given.

        If nothing is be found to match, return None.
        """
        return get_from_uuid(uuid)

    def get(self, uuid: str) -> Union[Molecule, Trajectory, Ensemble] | None:
        return self.entities.get(uuid)

    def prune(self) -> None:
        """
        Remove any entities that no longer exist in Blender
        """
        for uuid in list(self.entities):
            try:
                _ = self.entities[uuid].name
            except LinkedObjectError:
                self.remove_entity(uuid)

    @property
    def n_items(self) -> int:
        "The number of items being tracked by this session."
        return len(self.entities)

    def __repr__(self) -> str:
        return f"MNSession with {len(self.molecules)} molecules, {len(self.trajectories)} trajectories and {len(self.ensembles)} ensembles."

    def pickle(self, filepath) -> None:
        pickle_path = self.stashpath(filepath)

        make_paths_relative(self.trajectories)
        self.entities = trim(self.entities)

        # don't save anything if there is nothing to save
        if self.n_items == 0:
            return None

        with open(pickle_path, "wb") as f:
            pk.dump(self, f)

        print(f"Saved session to: {pickle_path}")

    def load(self, filepath) -> None:
        pickle_path = self.stashpath(filepath)
        if not os.path.exists(pickle_path):
            raise FileNotFoundError(f"MNSession file `{pickle_path}` not found")
        with open(pickle_path, "rb") as f:
            session = pk.load(f)

        # TODO: clear up in later versions
        # this handles reloading sessions which don't have the `entities` attribute
        if hasattr(session, "entities"):
            for item in session.entities.values():
                self.register_entity(item)
        else:
            for mol in session.molecules.values():
                self.register_entity(mol)

            for traj in session.trajectories.values():
                self.register_entity(traj)

            for ens in session.ensembles.values():
                self.register_entity(ens)

        print(f"Loaded a MNSession from: {pickle_path}")

    def stashpath(self, filepath) -> str:
        return f"{filepath}.MNSession"

    def clear(self) -> None:
        """Remove references to all molecules, trajectories and ensembles."""
        self.entities = {}

    def remove(self, uuid: str) -> None:
        """Remove an entity by uuid"""
        if uuid not in self.entities:
            raise ValueError(f"No entity with UUID '{uuid}'")
        entity = self.entities[uuid]
        bpy.data.objects.remove(entity.object, do_unlink=True)
        self.remove_entity(uuid)

    def add_trajectory(
        self,
        universe: mda.Universe | AtomGroup,
        style: str | None = "vdw",
        name: str = "NewUniverseObject",
    ) -> Trajectory:
        """
        Add a new trajectory

        Parameters
        ----------
        universe: mda.Universe | AtomGroup, required
            MDAnalysis Universe or AtomGroup instance

        style: str | None, optional
            The style to apply to the Universe or AtomGroup.

        name: str, optional
            Name of the trajectory object in Blender

        Returns
        -------
        Trajectory
            The newly added Trajectory instance

        """
        if style is not None and style not in styles_mapping:
            raise ValueError(
                f"Invalid style '{style}'. Supported styles are {[key for key in styles_mapping.keys()]}"
            )
        if isinstance(universe, AtomGroup):
            trajectory = Trajectory(universe.universe)
        else:
            trajectory = Trajectory(universe)
        trajectory.create_object(name=name, style=style)
        if isinstance(universe, AtomGroup):
            _selection = trajectory.add_selection_from_atomgroup(
                universe, name="InitialAtomGroupSelection"
            )
            # TODO: Add style only from the AtomGroup selection
            #       Will be fixed as part of style branch support for Trajectories
        return trajectory

    def get_trajectory(
        self,
        name: str,
    ) -> Trajectory:
        """
        Get trajectory instance by name

        Parameters
        ----------
        name: str, required
            Name of the trajectory object

        Returns
        -------
        Trajectory
            A Trajectory instance

        Raises
        ------
        ValueError if trajectory is not found

        """
        for v in self.entities.values():
            if isinstance(v, Trajectory) and v.object.name == name:
                return v
        raise ValueError(f"No trajectory named '{name}'")

    def remove_trajectory(
        self,
        trajectory: Trajectory | str,
    ) -> None:
        """
        Remove an existing trajectory

        Parameters
        ----------
        trajectory: Trajectory | str, required
            A Trajectory instance or name of the trajectory

        Returns
        -------
        None

        Raises
        ------
        ValueError if trajectory name is not found

        """
        instance = trajectory
        if isinstance(trajectory, str):
            instance = self.get_trajectory(trajectory)
        self.remove(instance.uuid)


def get_session(context: Context | None = None) -> MNSession:
    if not context:
        context = bpy.context
    return context.scene.MNSession


@persistent
def _pickle(filepath) -> None:
    get_session().pickle(filepath)


@persistent
def _load(filepath: str, printing: str = "quiet") -> None:
    """
    Load a session from the specified file path.

    This function attempts to load a session from the given file path using the
    `get_session().load(filepath)` method. If the file path is empty, the function
    returns immediately without attempting to load anything. If the file is not found,
    it handles the `FileNotFoundError` exception and optionally prints a message
    based on the `printing` parameter.

    Args:
        filepath (str): The path to the file from which to load the session. If this
            is an empty string, the function will return without doing anything.
        printing (str, optional): Controls the verbosity of the function. If set to
            "verbose", a message will be printed when the file is not found. Defaults
            to "quiet".

    Returns:
        None: This function does not return any value.

    Raises:
        FileNotFoundError: If the file specified by `filepath` does not exist and
            `printing` is set to "verbose", a message will be printed.
    """
    # the file hasn't been saved or we are opening a fresh file, so don't
    # attempt to load anything
    if filepath == "":
        return None
    try:
        get_session().load(filepath)
    except FileNotFoundError:
        if printing == "verbose":
            print("No MNSession found to load for this .blend file.")
        else:
            pass


class MN_OT_Session_Remove_Item(bpy.types.Operator):
    bl_idname = "mn.session_remove_item"
    bl_label = "Remove"
    bl_description = "Remove this item from the internal Molecular Nodes session"
    bl_options = {"REGISTER", "UNDO"}

    uuid: StringProperty()  # type: ignore

    def invoke(self, context: Context, event):
        session = get_session()

        return context.window_manager.invoke_confirm(
            self,
            event=event,
            title="Permanently delete item?",
            message=f"Any links to objects that rely upon this item will be lost.  {session.get(self.uuid)}",
        )

    def execute(self, context: Context):
        get_session().remove(self.uuid)

        return {"FINISHED"}


class MN_OT_Session_Create_Object(bpy.types.Operator):
    bl_idname = "mn.session_create_object"
    bl_label = "Create Object"
    bl_description = "Create a new object linked to this item"
    bl_options = {"REGISTER", "UNDO"}

    uuid: StringProperty()  # type: ignore

    def execute(self, context: Context):
        item = get_session().get(self.uuid)
        item.create_object()
        return {"FINISHED"}


CLASSES = [MN_OT_Session_Remove_Item, MN_OT_Session_Create_Object]
