import pickle as pk
import bpy
import os

from bpy.types import Context
from .io.parse.molecule import Molecule
from .io.md.universe import MNUniverse
from .io.parse.ensemble import Ensemble
from typing import List, Dict, Union
from bpy.app.handlers import persistent
from bpy.props import StringProperty, IntProperty, EnumProperty


class MN_OT_Session_Remove_Item(bpy.types.Operator):
    bl_idname = "mn.session_remove_item"
    bl_label = "Remove"
    bl_description = "Remove this item from the internal Molecular Nodes session"
    bl_options = {"REGISTER", "UNDO"}

    uuid: StringProperty()  # type: ignore

    def invoke(self, context: Context, event):
        session = context.scene.MNSession

        return context.window_manager.invoke_confirm(
            self,
            event=event,
            title="Permanently delete item?",
            message=f"Any links to objects that rely upon this item will be lost.  {session.get(self.uuid)}",
        )

    def execute(self, context: Context):
        session = context.scene.MNSession
        session.remove(self.uuid)

        return {"FINISHED"}


class MN_OT_Session_Create_Model(bpy.types.Operator):
    bl_idname = "mn.session_create_model"
    bl_label = "Create Model"
    bl_description = "Create a new model object linked to this item"
    bl_options = {"REGISTER", "UNDO"}

    uuid: StringProperty()  # type: ignore

    def execute(self, context: Context):
        session = context.scene.MNSession
        item = session.get(self.uuid)
        item.create_model()
        return {"FINISHED"}


def trim(dictionary: dict):
    to_pop = []
    for name, item in dictionary.items():
        try:
            if isinstance(item.object, bpy.types.Object):
                item.name = item.object.name
                item.object = None
            if hasattr(item, "frames"):
                if isinstance(item.frames, bpy.types.Collection):
                    item.frames_name = item.frames.name
                    item.frames = None

        except ReferenceError as e:
            to_pop.append(name)
            print(
                Warning(
                    f"Object reference for {item} broken, removing this item from the session: `{e}`"
                )
            )

    for name in to_pop:
        dictionary.pop(name)
    return dictionary


class MNSession:
    def __init__(self) -> None:
        self.molecules: Dict[str, Molecule] = {}
        self.universes: Dict[str, MNUniverse] = {}
        self.ensembles: Dict[str, Ensemble] = {}

    def items(self):
        "Return UUID and item for all molecules, universes and ensembles being tracked."
        return (
            list(self.molecules.items())
            + list(self.universes.items())
            + list(self.ensembles.items())
        )

    def get_object(self, uuid: str) -> bpy.types.Object | None:
        """
        Try and get an object from Blender's object database that matches the uuid given.

        If nothing is be found to match, return None.
        """
        for bob in bpy.data.objects:
            try:
                if bob.mn.uuid == uuid:
                    return bob
            except Exception as e:
                print(e)

        return None

    def remove(self, uuid: str) -> None:
        "Remove the item from the list."
        self.molecules.pop(uuid, None)
        self.universes.pop(uuid, None)
        self.ensembles.pop(uuid, None)

    def get(self, uuid: str) -> Union[Molecule, MNUniverse, Ensemble]:
        for id, item in self.items():
            if item.uuid == uuid:
                return item

        return None

    @property
    def n_items(self) -> int:
        "The number of items being tracked by this session."
        length = 0

        for dic in [self.molecules, self.universes, self.ensembles]:
            length += len(dic)
        return length

    def __repr__(self) -> str:
        return f"MNSession with {len(self.molecules)} molecules, {len(self.universes)} universes and {len(self.ensembles)} ensembles."

    def pickle(self, filepath) -> None:
        pickle_path = self.stashpath(filepath)

        self.molecules = trim(self.molecules)
        self.universes = trim(self.universes)
        self.ensembles = trim(self.ensembles)

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

        for uuid, item in session.items():
            item.object = bpy.data.objects[item.name]
            if hasattr(item, "frames") and hasattr(item, "frames_name"):
                item.frames = bpy.data.collections[item.frames_name]

        for uuid, mol in session.molecules.items():
            self.molecules[uuid] = mol

        for uuid, uni in session.universes.items():
            self.universes[uuid] = uni

        for uuid, ens in session.ensembles.items():
            self.ensembles[uuid] = ens

        print(f"Loaded a MNSession from: {pickle_path}")

    def stashpath(self, filepath) -> str:
        return f"{filepath}.MNSession"

    def clear(self) -> None:
        """Remove references to all molecules, universes and ensembles."""
        # for mol in self.molecules:
        #     try:
        #         o = mol.object
        #         mol.object = None
        #         bpy.data.objects.remove(o)
        #         if mol.frames is not None:
        #             for obj in mol.frames.objects:
        #                 bpy.data.objects.remove(obj)
        #             c = mol.frames
        #             mol.frames = None
        #             bpy.data.collections.remove(c)
        #     except ReferenceError:
        #         pass
        # for univ in self.universes:
        #     try:
        #         bpy.data.objects.remove(univ.object)
        #     except ReferenceError:
        #         pass
        # for ens in self.ensembles:
        #     try:
        #         bpy.data.objects.remove(ens.object)
        #     except ReferenceError:
        #         pass

        self.molecules.clear()
        self.universes.clear()
        self.ensembles.clear()


@persistent
def _pickle(filepath) -> None:
    bpy.context.scene.MNSession.pickle(filepath)


@persistent
def _load(filepath) -> None:
    # the file hasn't been saved or we are opening a fresh file, so don't
    # attempt to load anything
    if filepath == "":
        return None
    try:
        bpy.context.scene.MNSession.load(filepath)
    except FileNotFoundError:
        print("No MNSession found to load for this .blend file.")


CLASSES = [MN_OT_Session_Remove_Item, MN_OT_Session_Create_Model]
