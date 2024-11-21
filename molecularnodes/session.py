import os
import pickle as pk
from typing import Dict, Union

import bpy
from pathlib import Path
from bpy.app.handlers import persistent
from bpy.props import StringProperty
from bpy.types import Context

from .entities.ensemble.ensemble import Ensemble
from .entities.molecule.molecule import Molecule
from .entities.trajectory.trajectory import Trajectory


def path_relative_to_blend_wd(filepath: str | Path) -> Path:
    "Get the path of something, relative to the working directory of the current .blend file"
    blend_working_directory = bpy.path.abspath("//")
    if blend_working_directory == "":
        raise ValueError(
            "Unable to get current working directly, .blend file not saved"
        )

    return Path(filepath).relative_to(Path(blend_working_directory))


def make_paths_relative(trajectories: Dict[str, Trajectory]) -> None:
    for key, traj in trajectories.items():
        traj.universe.load_new(
            path_relative_to_blend_wd(traj.universe.trajectory.filename)
        )
        traj.save_filepaths_on_object()


def find_matching_object(uuid):
    for obj in bpy.data.objects:
        if obj.mn.uuid == uuid:
            return obj

    return None


class MNSession:
    def __init__(self) -> None:
        self.entities: Dict[str, Molecule | Trajectory | Ensemble] = {}

    @property
    def molecules(self) -> dict:
        return {
            key: mol for key, mol in self.entities.items() if isinstance(mol, Molecule)
        }

    @property
    def trajectories(self) -> dict:
        return {
            key: traj
            for key, traj in self.entities.items()
            if isinstance(traj, Trajectory)
        }

    @property
    def ensembles(self) -> dict:
        return {
            key: ens for key, ens in self.entities.items() if isinstance(ens, Ensemble)
        }

    def get(self, uuid: str) -> Union[Molecule, Trajectory, Ensemble]:
        return self.entities.get(uuid)

    def __repr__(self) -> str:
        return f"MNSession with {len(self.molecules)} molecules, {len(self.trajectories)} trajectories and {len(self.ensembles)} ensembles."

    def __len__(self) -> int:
        return len(self.entities)

    def trim(self) -> None:
        to_pop = []
        for name, item in self.entities.items():
            # currently there are problems with pickling the functions so we have to just
            # clean up any calculations that are created on saving. Could potentially convert
            # it to a string and back but that is likely a job for better implementations
            if hasattr(item, "calculations"):
                item.calculations = {}

            if item.object is None:
                to_pop.append(name)

        for name in to_pop:
            self.entities.pop(name)

    def pickle(self, filepath) -> None:
        path = Path(filepath)
        self.trim()
        if len(self) == 0:
            return None

        make_paths_relative(self.trajectories)

        # don't save anything if there is nothing to save
        if len(self) == 0:
            # if we aren't saving anything, remove the currently existing session file
            # so that it isn't reloaded when we load the save with old session information
            if path.exists() and path.suffix == ".MNSession":
                os.remove(filepath)
            return None

        with open(filepath, "wb") as f:
            pk.dump(self, f)

        print(f"Saved MNSession to: {filepath}")

    def load(self, filepath) -> None:
        "Load all of the entities from a previously saved MNSession"
        path = Path(filepath)

        if not path.exists():
            raise FileNotFoundError(f"MNSession file `{path}` not found")

        with open(path, "rb") as f:
            loaded_session = pk.load(f)
        current_session = bpy.context.scene.MNSession

        # merge the loaded session with current session, handling if they used the old
        # structure of separating entities into different categories
        if hasattr(loaded_session, "entities"):
            current_session.entites + loaded_session.entities
        else:
            items = []
            for attr in ["molecules", "trajectories", "ensembles"]:
                try:
                    items.append((getattr(loaded_session, attr)))
                except AttributeError:
                    pass
            print(f"{items=}")
            for key, item in items:
                current_session.entities[key] = item

        print(f"Loaded a MNSession from: {filepath}")

    def stashpath(self, filepath) -> str:
        return f"{filepath}.MNSession"

    def clear(self) -> None:
        """Remove references to all molecules, trajectories and ensembles."""
        self.entities.clear()


def get_session(context: Context | None = None) -> MNSession:
    if isinstance(context, Context):
        return context.scene.MNSession
    else:
        return bpy.context.scene.MNSession


@persistent
def _pickle(filepath) -> None:
    session = get_session()
    session.pickle(session.stashpath(filepath))


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
        session = get_session()
        session.load(session.stashpath(filepath))
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
