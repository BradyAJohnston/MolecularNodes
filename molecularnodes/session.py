import pickle as pk
import bpy
import os
from .io.parse.molecule import Molecule
from .io.parse.mda import MNUniverse
from .io.parse.ensemble import Ensemble
from typing import List, Dict
from bpy.app.handlers import persistent


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
            print(f"Reference to {item} broken, not saving. {e}")

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

    def __repr__(self) -> str:
        return f"MNSession with {len(self.molecules)} molecules, {len(self.universes)} universes and {len(self.ensembles)} ensembles."

    def pickle(self, filepath) -> None:
        pickle_path = self.stashpath(filepath)

        self.molecules = trim(self.molecules)
        self.universes = trim(self.universes)
        self.ensembles = trim(self.ensembles)

        if len(self.universes) == 0:
            print(f"Skipping saving of MNSession, {len(self.universes)=}")
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
            if hasattr(item, "frames"):
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
