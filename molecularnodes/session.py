import pickle as pk
import warnings
from contextlib import chdir
from pathlib import Path
from typing import Dict, Union
import bpy
from bpy.app.handlers import persistent
from bpy.props import StringProperty
from bpy.types import Context
from databpy.object import LinkedObjectError, get_from_uuid
from .entities import Molecule
from .entities.base import MolecularEntity
from .entities.ensemble.base import Ensemble
from .entities.reload import can_reload, reload_entity


def trim(dictionary: dict):
    """Drop entities whose Blender object no longer exists.

    ``.object`` raises ``LinkedObjectError`` rather than returning ``None`` once
    the object is gone, so this has to catch as well as test - otherwise a single
    entity deleted from the outliner takes the whole session save down with it,
    and every other entity in the file silently loses its state.
    """
    dic = dictionary.copy()
    for key in list(dic.keys()):
        try:
            if dic[key].object is None:
                dic.pop(key)
        except LinkedObjectError:
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


def _remap_trajectory_paths(trajectories: Dict[str, Molecule], remap) -> None:
    for key, traj in trajectories.items():
        if not _has_ondisk_trajectory(traj):
            continue
        # multi-file universes (ChainReader) list every coordinate file in
        # `filenames`, while `filename` only holds the first one
        filenames = getattr(traj.universe.trajectory, "filenames", None)
        # save linked universe frame
        uframe = traj.uframe
        if filenames is not None:
            traj.universe.load_new([remap(f) for f in filenames])
        else:
            traj.universe.load_new(remap(traj.universe.trajectory.filename))
        # restore linked universe frame
        traj.uframe = uframe
        traj._save_filepaths_on_object()


def _make_trajectory_paths_relative(trajectories: Dict[str, Molecule]) -> None:
    _remap_trajectory_paths(trajectories, _make_path_relative)


def _make_trajectory_paths_absolute(trajectories: Dict[str, Molecule]) -> None:
    _remap_trajectory_paths(trajectories, _make_path_absolute)


def _make_path_relative(filepath):
    "Take a path and make it relative"
    try:
        return Path(filepath).relative_to(Path.cwd(), walk_up=True)
    except ValueError:
        return filepath


def _make_path_absolute(filepath):
    "Take a path and make it absolute"
    try:
        return Path(filepath).resolve()
    except ValueError:
        return filepath


class MNSession:
    def __init__(self) -> None:
        self.entities: Dict[str, Union[Molecule, Ensemble]] = {}

    @property
    def molecules(self) -> Dict[str, Molecule]:
        return {k: v for k, v in self.entities.items() if isinstance(v, Molecule)}

    @property
    def ensembles(self) -> Dict[str, Ensemble]:
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
        return f"MNSession with {len(self.molecules)} molecules and {len(self.ensembles)} ensembles."

    def pickle(self, filepath) -> None:
        pickle_path = self.stashpath(filepath)

        self.entities = trim(self.entities)

        # don't save anything if there is nothing to save
        if self.n_items == 0:
            return None

        # skip entities that can't be pickled, so that one bad entity doesn't lose
        # the session state for everything else in the scene
        picklable = {}
        for uuid, entity in self.entities.items():
            try:
                pk.dumps(entity)
                picklable[uuid] = entity
            except Exception as e:
                warnings.warn(
                    f"Not saving `{entity.name}` with the session, "
                    f"as it cannot be serialized: {e}"
                )
        if len(picklable) == 0:
            return None

        all_entities = self.entities
        self.entities = picklable
        _make_trajectory_paths_relative(self.molecules)
        # dump to a temporary file and swap it into place, so that a failed dump
        # can't leave a truncated session file next to the .blend
        temp_path = pickle_path.with_name(f"{pickle_path.name}.tmp")
        try:
            with open(temp_path, "wb") as f:
                pk.dump(self, f)
            temp_path.replace(pickle_path)
        finally:
            temp_path.unlink(missing_ok=True)
            _make_trajectory_paths_absolute(self.molecules)
            self.entities = all_entities

        print(f"Saved session to: {pickle_path}")

    def load(self, filepath) -> None:
        pickle_path = self.stashpath(filepath)
        if not pickle_path.exists():
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

        _make_trajectory_paths_absolute(self.molecules)

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

    def stashpath(self, filepath) -> Path:
        return Path(f"{filepath}.MNSession")

    def clear(self) -> None:
        """Remove references to all molecules, trajectories and ensembles."""
        self.entities = {}

    def remove(self, uuid: str) -> None:
        """
        Remove an entity by uuid, along with the objects backing it.

        An entity whose object has already been deleted - from the outliner, or
        with `bpy.data.objects.remove()` - is still dropped from the session
        rather than raising, so that removing is always safe to call.
        """
        if uuid not in self.entities:
            raise ValueError(f"No entity with UUID '{uuid}'")
        entity = self.entities[uuid]
        obj = self.get_object(uuid)
        if obj is not None:
            bpy.data.objects.remove(obj, do_unlink=True)
        if hasattr(entity, "annotations"):
            bob = entity.annotations.bob
            if bob is not None and self.get_object(bob.uuid) is not None:
                bpy.data.objects.remove(bob.object, do_unlink=True)
        self.remove_entity(uuid)


def get_session(context: Context | None = None) -> MNSession:
    if not context:
        context = bpy.context
    return context.scene.MNSession


def get_entity(context: Context | None = None) -> Molecule | Ensemble:
    session = get_session(context)
    return session.match(context.active_object)


@persistent
def _pickle(filepath) -> None:
    with chdir(Path(filepath).parent):
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
    # the objects belonging to the file we just closed are gone, so drop the
    # entities that pointed at them before registering this file's own -
    # otherwise they linger in the session and raise LinkedObjectError
    session = get_session()
    session.prune()

    # the file hasn't been saved or we are opening a fresh file, so don't
    # attempt to load anything
    if filepath == "":
        return None
    try:
        with chdir(Path(filepath).parent):
            session.load(filepath)
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
