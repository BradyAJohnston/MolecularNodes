import os
import pickle as pk
from typing import Dict, Union

import bpy
from pathlib import Path
from bpy.app.handlers import persistent
from bpy.props import StringProperty
from bpy.types import Context, Operator

from .entities import Ensemble, Molecule, Trajectory


def path_relative_to_blend(target_path: str | Path) -> Path:
    """Get a path relative to the current .blend file"""
    blend_path = bpy.data.filepath
    if blend_path == "":
        raise ValueError(
            ".blend file has not yet been saved, unable to get relative path"
        )

    blender_folder = Path(blend_path).parent.absolute()

    target_path = Path(target_path)
    if not target_path.is_absolute():
        target_path = (blender_folder / target_path).resolve()

    # Get the relative path
    try:
        relative_path = Path(os.path.relpath(target_path, blender_folder))
        return relative_path
    except ValueError as e:
        # Handle case where paths are on different drives (Windows)
        return target_path


def make_paths_relative(trajectories: Dict[str, Trajectory]) -> None:
    for key, traj in trajectories.items():
        newpath = path_relative_to_blend(traj.universe.trajectory.filename)
        cwd = Path.cwd()
        try:
            os.chdir(Path(bpy.data.filepath).parent)
            traj.universe.load_new(newpath)
            traj.save_filepaths_on_object()
        finally:
            os.chdir(cwd)


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
            loaded_session: MNSession = pk.load(f)
            if not isinstance(loaded_session, MNSession):
                raise ValueError(
                    f"Loaded .pkl object is not a MNSession, instead: {loaded_session=}"
                )

        current_session: MNSession = bpy.context.scene.MNSession

        # merge the loaded session with current session, handling if they used the old
        # structure of separating entities into different categories
        if hasattr(loaded_session, "entities"):
            current_session.entities | loaded_session.entities
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


class MN_OT_Session_Remove_Item(Operator):
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


class MN_OT_Session_Create_Object(Operator):
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
