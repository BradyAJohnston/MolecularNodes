import pickle as pk
import bpy
import os
from .io.parse.molecule import Molecule
from .io.parse.mda import MNUniverse
from .io.parse.ensemble import Ensemble
from typing import List
from bpy.app.handlers import persistent


def trim(items: list):
    trimmed_items = []
    for item in items:
        try:
            if isinstance(item.object, bpy.types.Object):
                item.name = item.object.name
                item.object = None
            if hasattr(item, "frames"):
                if isinstance(item.frames, bpy.types.Collection):
                    item.frames_name = item.frames.name
                    item.frames = None
            trimmed_items.append(item)

        except ReferenceError as e:
            print(f"Reference to {item} broken, not saving. {e}")
    return trimmed_items


class MNSession:
    def __init__(self) -> None:
        self.molecules: List[Molecule] = []
        self.universes: List[MNUniverse] = []
        self.ensembles: List[Ensemble] = []

    def lists(self):
        return [self.molecules, self.universes, self.ensembles]

    def __repr__(self) -> str:
        return f"MNSession with {len(self.molecules)} molecules, {len(self.universes)} universes and {len(self.ensembles)} ensembles."

    def pickle(self, filepath) -> None:
        pickle_path = self.stashpath(filepath)

        self.molecules = trim(self.molecules)
        self.universes = trim(self.universes)
        self.ensembles = trim(self.ensembles)

        with open(pickle_path, "wb") as f:
            pk.dump(self, f)

        print(f"Saved session to: {pickle_path}")

    def load(self, filepath) -> None:
        pickle_path = self.stashpath(filepath)
        if not os.path.exists(pickle_path):
            raise FileNotFoundError(f"MNSession file `{pickle_path}` not found")
        with open(pickle_path, "rb") as f:
            session = pk.load(f)

        for items in session.lists():
            for item in items:
                item.object = bpy.data.objects[item.name]
                if hasattr(item, "frames"):
                    item.frames = bpy.data.collections[item.frames_name]

        for item in session.molecules:
            self.molecules.append(item)

        for item in session.universes:
            self.universes.append(item)

        for item in session.ensembles:
            self.ensembles.append(item)

        print(f"Loaded a MNSession from: {pickle_path}")

    def stashpath(self, filepath) -> str:
        return f"{filepath}.MNSession"

    def clean(self) -> None:
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

        self.molecules = []
        self.universes = []
        self.ensembles = []


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
