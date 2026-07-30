import os
import pickle as pk
from contextlib import chdir
from pathlib import Path
from typing import Dict, Union
import bpy
import MDAnalysis as mda
from bpy.app.handlers import persistent
from bpy.props import StringProperty
from bpy.types import Context
from databpy.object import LinkedObjectError, get_from_uuid
from MDAnalysis.core.groups import AtomGroup
from .entities import Molecule
from .entities.base import MolecularEntity
from .entities.ensemble.base import Ensemble
from .entities.reload import can_reload, reload_entity
from .nodes.nodes import styles_mapping


def trim(dictionary: dict):
    dic = dictionary.copy()
    for key in list(dic.keys()):
        if dic[key].object is None:
            dic.pop(key)
    return dic


def _has_ondisk_trajectory(traj: Molecule) -> bool:
    """Whether a trajectory is backed by a relocatable on-disk trajectory file.

    Streaming (IMD) trajectories and in-memory universes (e.g. biotite-converted
    molecules, whose ``filename`` is a ``BiotiteWrapper`` or ``None``) have no
    on-disk trajectory file whose path can be made relative/absolute.
    """
    filename = traj.universe.trajectory.filename
    if not isinstance(filename, (str, Path)):
        return False
    return not str(filename).startswith("imd://")


def _make_trajectory_paths_relative(trajectories: Dict[str, Molecule]) -> None:
    for key, traj in trajectories.items():
        if not _has_ondisk_trajectory(traj):
            continue
        # save linked universe frame
        uframe = traj.uframe
        traj.universe.load_new(_make_path_relative(traj.universe.trajectory.filename))
        # restore linked universe frame
        traj.uframe = uframe
        traj._save_filepaths_on_object()


def _make_trajectory_paths_absolute(trajectories: Dict[str, Molecule]) -> None:
    for key, traj in trajectories.items():
        if not _has_ondisk_trajectory(traj):
            continue
        # save linked universe frame
        uframe = traj.uframe
        traj.universe.load_new(_make_path_absolute(traj.universe.trajectory.filename))
        # restore linked universe frame
        traj.uframe = uframe
        traj._save_filepaths_on_object()


def _make_path_relative(filepath):
    "Take a path and make it relative"
    try:
        return os.path.relpath(filepath)
    except ValueError:
        return filepath


def _make_path_absolute(filepath):
    "Take a path and make it absolute"
    try:
        return os.path.abspath(filepath)
    except ValueError:
        return filepath


class MNSession:
    def __init__(self) -> None:
        self.entities: Dict[str, Union[Molecule, Ensemble]] = {}

    @property
    def molecules(self) -> Dict[str, Molecule]:
        return {k: v for k, v in self.entities.items() if isinstance(v, Molecule)}

    @property
    def trajectories(self) -> Dict[str, Molecule]:
        # return a filtered dictionary of only the trajectories using isinstance(item, Molecule)
        return {k: v for k, v in self.entities.items() if isinstance(v, Molecule)}

    @property
    def ensembles(self) -> Dict[str, Ensemble]:
        # return a filtered dictionary of only the ensembles using isinstance(item, Ensemble)
        return {k: v for k, v in self.entities.items() if isinstance(v, Ensemble)}

    def register_entity(self, item: MolecularEntity) -> None:
        """Add entity to the session, keyed by its uuid."""
        self.entities[item.uuid] = item

    def remove_entity(self, uuid: str) -> None:
        """Remove entity from the session."""
        del self.entities[uuid]

    def match(self, obj: bpy.types.Object) -> Union[Molecule, Ensemble]:
        return self.get(obj.uuid)

    def get_object(self, uuid: str) -> bpy.types.Object | None:
        """
        Try and get an object from Blender's object database that matches the uuid given.

        If nothing is be found to match, return None.
        """
        try:
            return get_from_uuid(uuid)
        except LinkedObjectError:
            return None

    def get(self, uuid: str) -> Union[Molecule, Ensemble] | None:
        return self.entities.get(uuid)

    def prune(self) -> None:
        """
        Remove any session entities whose Blender object no longer exists.
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

        self.entities = trim(self.entities)

        # don't save anything if there is nothing to save
        if self.n_items == 0:
            return None

        _make_trajectory_paths_relative(self.trajectories)

        with open(pickle_path, "wb") as f:
            pk.dump(self, f)

        _make_trajectory_paths_absolute(self.trajectories)

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

        _make_trajectory_paths_absolute(self.trajectories)

        print(f"Loaded a MNSession from: {pickle_path}")

    def remove_draw_handlers(self) -> None:
        for e in self.entities.values():
            if hasattr(e, "annotations"):
                e.annotations._draw_handler_remove()

    def add_draw_handlers(self) -> None:
        self.prune()
        for e in self.entities.values():
            if hasattr(e, "annotations") and e.annotations.visible:
                e.annotations._draw_handler_add()

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
        if hasattr(entity, "annotations"):
            if entity.annotations.bob:
                bpy.data.objects.remove(entity.annotations.bob.object, do_unlink=True)
        self.remove_entity(uuid)

    def add_trajectory(
        self,
        universe: mda.Universe | AtomGroup,
        style: str | None = "spheres",
        name: str = "NewUniverseObject",
    ) -> Molecule:
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
        Molecule
            The newly added Molecule instance

        """
        if style is not None and style not in styles_mapping:
            raise ValueError(
                f"Invalid style '{style}'. Supported styles are {[key for key in styles_mapping.keys()]}"
            )
        selection = None
        if isinstance(universe, AtomGroup):
            traj = Molecule(universe.universe, name=name)  # AtomGroup universe
            selection = universe  # AtomGroup
        else:
            traj = Molecule(universe, name=name)  # Universe
        traj.add_style(style=style, selection=selection)
        return traj

    def get_trajectory(
        self,
        name: str,
    ) -> Molecule:
        """
        Get trajectory instance by name

        Parameters
        ----------
        name: str, required
            Name of the trajectory object

        Returns
        -------
        Molecule
            A Molecule instance

        Raises
        ------
        ValueError if trajectory is not found

        """
        for v in self.entities.values():
            if isinstance(v, Molecule) and v.object.name == name:
                return v
        raise ValueError(f"No trajectory named '{name}'")

    def remove_trajectory(
        self,
        trajectory: Molecule | str,
    ) -> None:
        """
        Remove an existing trajectory

        Parameters
        ----------
        trajectory: Molecule | str, required
            A Molecule instance or name of the trajectory

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


def get_entity(context: Context | None = None) -> Molecule | Ensemble:
    session = get_session(context)
    return session.match(context.active_object)


@persistent
def _pickle(filepath) -> None:
    with chdir(os.path.dirname(filepath)):
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
        with chdir(os.path.dirname(filepath)):
            get_session().load(filepath)
    except FileNotFoundError:
        if printing == "verbose":
            print("No MNSession found to load for this .blend file.")
        else:
            pass


@persistent
def _remove_draw_handlers(filepath: str) -> None:
    get_session().remove_draw_handlers()


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


class MN_OT_Session_Reload_Item(bpy.types.Operator):
    bl_idname = "mn.session_reload_item"
    bl_label = "Reload"
    bl_description = (
        "Reload this entity's data into the session from its source file or PDB "
        "code, relinking the object to a live Molecular Nodes entity"
    )
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context: Context) -> bool:
        obj = context.active_object
        return obj is not None and can_reload(obj)

    def execute(self, context: Context):
        obj = context.active_object
        try:
            reload_entity(obj)
        except Exception as e:
            self.report({"ERROR"}, f"Failed to reload entity: {e}")
            return {"CANCELLED"}
        return {"FINISHED"}


CLASSES = [
    MN_OT_Session_Remove_Item,
    MN_OT_Session_Reload_Item,
]
