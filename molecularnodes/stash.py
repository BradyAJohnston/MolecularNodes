import pickle
import bpy
from .io.parse.molecule import Molecule
from .io.parse.mda import MNUniverse
from .io.parse.ensemble import Ensemble
from typing import List
from bpy.app.handlers import persistent
from pathlib import Path


class MNSession:
    def __init__(self) -> None:
        self.molecules: List[Molecule] = []
        self.universes: List[MNUniverse] = []
        self.ensembles: List[Ensemble] = []

    def lists(self):
        return [self.molecules, self.universes, self.ensembles]

    def load(self, session):
        for item in session.molecules:
            self.molecules.append(item)

        for item in session.universes:
            self.universes.append(item)

        for item in session.ensembles:
            self.ensembles.append(item)

    def __repr__(self) -> str:
        return f"MNSession with {self.molecules=}{self.universes=}{self.ensembles=}"


def stashpath():
    return f"{bpy.data.filepath}.MNSession"


@persistent
def _stash_save(scene) -> None:
    path = stashpath()
    stash_database(session=bpy.context.scene.MNSession, filepath=stashpath())
    print(f"Save database to: {path}")


@persistent
def _stash_load(scene) -> None:
    apply_database(stashpath())


def stash_database(session, filepath) -> None:
    # have to unlink
    for items in session.lists():
        for item in items:
            if isinstance(item.object, bpy.types.Object):
                if item.object is not None:
                    item.name = item.object.name
                    item.object = None

    with open(stashpath(), "wb") as f:
        pickle.dump(session, f)


def apply_database(file: Path) -> None:
    with open(file, "rb") as f:
        session = pickle.load(f)

    for items in session.lists():
        for item in items:
            item.object = bpy.data.objects[item.name]

    bpy.context.scene.MNSession.load(session)

    # return database
