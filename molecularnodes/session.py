import pickle as pk
import bpy
from .io.parse.molecule import Molecule
from .io.parse.mda import MNUniverse
from .io.parse.ensemble import Ensemble
from typing import List
from bpy.app.handlers import persistent


class MNSession:
    def __init__(self) -> None:
        self.molecules: List[Molecule] = []
        self.universes: List[MNUniverse] = []
        self.ensembles: List[Ensemble] = []

    def lists(self):
        return [self.molecules, self.universes, self.ensembles]

    def __repr__(self) -> str:
        return f"MNSession with {len(self.molecules)} molecules, {len(self.universes)} universes and {len(self.ensembles)} ensembles."

    def pickle(self) -> None:
        # have to unlink
        for items in self.lists():
            for item in items:
                if isinstance(item.object, bpy.types.Object):
                    if item.object is not None:
                        item.name = item.object.name
                        item.object = None

        with open(self.stashpath(), "wb") as f:
            pk.dump(self, f)

        print(f"Saved session to: {self.stashpath()}")

    def load(self) -> None:
        with open(self.stashpath(), "rb") as f:
            session = pk.load(f)

        for items in session.lists():
            for item in items:
                item.object = bpy.data.objects[item.name]

        for item in session.molecules:
            self.molecules.append(item)

        for item in session.universes:
            self.universes.append(item)

        for item in session.ensembles:
            self.ensembles.append(item)

    def stashpath(self) -> str:
        return f"{bpy.data.filepath}.MNSession"


@persistent
def _pickle(filepath) -> None:
    bpy.context.scene.MNSession.pickle()


@persistent
def _load(filepath):
    bpy.context.scene.MNSession.load()
