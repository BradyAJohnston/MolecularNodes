import os
import warnings
from pathlib import Path
from typing import Dict, Union
import bpy
from bpy.app.handlers import persistent
from bpy.props import StringProperty
from bpy.types import Context
from databpy.object import LinkedObjectError, get_from_uuid
from .blender.registered import object_props, object_uuid, scene_session
from .entities import Molecule
from .entities.density.base import Density
from .entities.ensemble.base import Ensemble
from .entities.reload import can_reload, reload_entity, restore_entity

# the source paths recorded on an entity object, which saving makes relative to
# the .blend so that a project folder can be moved or shared as a whole
_PATH_PROPS = ("filepath", "filepath_topology", "filepath_trajectory")


def _relative_to(path: str, folder: Path) -> str:
    """``path`` as a Blender-relative ``//`` path from ``folder``.

    The path is returned unchanged when it is empty, or when no relative path to
    it exists (a different drive on Windows).
    """
    if not path:
        return path
    try:
        relative = os.path.relpath(bpy.path.abspath(path), folder)
    except ValueError:
        return path
    return "//" + Path(relative).as_posix()


def _make_paths_relative(folder: Path) -> None:
    """Make the source paths on every entity object relative to ``folder``."""
    for obj in bpy.data.objects:
        # linked objects can't be edited, and their paths are relative to
        # their own library file
        if obj.library is not None:
            continue
        mn = object_props(obj)
        if mn.entity_type == "None":
            continue
        for name in _PATH_PROPS:
            setattr(mn, name, _relative_to(getattr(mn, name), folder))
        filepaths = mn.filepaths_trajectory
        if filepaths:
            mn.filepaths_trajectory = [_relative_to(f, folder) for f in filepaths]


def _restore_entities() -> None:
    """Rebuild the session entity of every entity object in the file.

    Everything needed to rebuild an entity is recorded on its object, so each one
    is reloaded from its source files (or PDB code) and linked back to the object
    under the uuid it already has. One entity failing to reload doesn't stop the
    rest.
    """
    for obj in bpy.data.objects:
        if obj.library is not None or not can_reload(obj):
            continue
        try:
            restore_entity(obj)
        except Exception as e:
            warnings.warn(
                f"Could not restore `{obj.name}` from its source files, so it "
                "won't update with the scene. Once its source files are "
                "available, select it and use Reload in the Molecular Nodes "
                f"panel to relink it. {e}"
            )


class MNSession:
    def __init__(self) -> None:
        self.entities: Dict[str, Union[Molecule, Ensemble, Density]] = {}

    @property
    def molecules(self) -> Dict[str, Molecule]:
        return {k: v for k, v in self.entities.items() if isinstance(v, Molecule)}

    @property
    def ensembles(self) -> Dict[str, Ensemble]:
        return {k: v for k, v in self.entities.items() if isinstance(v, Ensemble)}

    def register_entity(self, item: Molecule | Ensemble | Density) -> None:
        """Add entity to the session, keyed by its uuid."""
        self.entities[item.uuid] = item

    def remove_entity(self, uuid: str) -> None:
        """Remove entity from the session."""
        del self.entities[uuid]

    def match(self, obj: bpy.types.Object) -> Molecule | Ensemble | Density | None:
        return self.get(object_uuid(obj))

    def get_object(self, uuid: str) -> bpy.types.Object | None:
        """
        Try and get an object from Blender's object database that matches the uuid given.

        If nothing is be found to match, return None.
        """
        try:
            return get_from_uuid(uuid)
        except LinkedObjectError:
            return None

    def get(self, uuid: str) -> Union[Molecule, Ensemble, Density] | None:
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

    def remove_draw_handlers(self) -> None:
        for e in self.entities.values():
            annotations = getattr(e, "annotations", None)
            if annotations is not None:
                annotations._draw_handler_remove()

    def add_draw_handlers(self) -> None:
        self.prune()
        for e in self.entities.values():
            annotations = getattr(e, "annotations", None)
            if annotations is not None and annotations.visible:
                annotations._draw_handler_add()

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
        annotations = getattr(entity, "annotations", None)
        if annotations is not None:
            bob = annotations.bob
            if bob is not None and self.get_object(bob.uuid) is not None:
                bpy.data.objects.remove(bob.object, do_unlink=True)
        self.remove_entity(uuid)


def get_session(context: Context | None = None) -> MNSession:
    if not context:
        context = bpy.context
    assert context.scene is not None
    return scene_session(context.scene)


def get_entity(context: Context | None = None) -> Molecule | Ensemble:
    if not context:
        context = bpy.context
    session = get_session(context)
    assert context.active_object is not None
    # the active object can also be a Density or have no entity; the narrower
    # return type is kept for the UI callers
    return session.match(context.active_object)  # ty: ignore[invalid-return-type]


@persistent
def _save_pre(filepath: str) -> None:
    """Make entity source paths relative to the .blend about to be written."""
    _make_paths_relative(Path(filepath).parent)


@persistent
def _load(filepath: str) -> None:
    """Rebuild the session entities of the file that was just opened.

    The entities are rebuilt from what is recorded on their objects. A
    ``.MNSession`` file saved next to the .blend by older versions is never read:
    it is a pickle, and unpickling a file that came with a downloaded .blend
    could run arbitrary code.
    """
    # the objects belonging to the file we just closed are gone, so drop the
    # entities that pointed at them before registering this file's own -
    # otherwise they linger in the session and raise LinkedObjectError
    session = get_session()
    session.prune()

    # the file hasn't been saved or we are opening a fresh file, so there is
    # nothing to restore
    if filepath == "":
        return None
    _restore_entities()
    if Path(f"{filepath}.MNSession").exists():
        print(
            f"Ignoring `{Path(filepath).name}.MNSession`: Molecular Nodes no longer "
            "reads these files, and rebuilds entities from their source files "
            "instead. It can be deleted."
        )


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

        assert context.window_manager is not None
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
        assert obj is not None
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
