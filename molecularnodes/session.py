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
        return f"MNSession with {len(self.molecules)} molecules, {len(self.universes)} universes and {len(self.ensembles)} ensembles."


def stashpath():
    return f"{bpy.data.filepath}.MNSession"


@persistent
def _session_pickle(scene) -> None:
    path = stashpath()
    session_pickle(session=bpy.context.scene.MNSession, filepath=stashpath())
    print(f"Save database to: {path}")


@persistent
def _session_load(scene) -> None:
    session_load(stashpath())


def session_pickle(session, filepath) -> None:
    # have to unlink
    for items in session.lists():
        for item in items:
            if isinstance(item.object, bpy.types.Object):
                if item.object is not None:
                    item.name = item.object.name
                    item.object = None

    with open(stashpath(), "wb") as f:
        pickle.dump(session, f)


def session_load(file: Path) -> None:
    with open(file, "rb") as f:
        session = pickle.load(f)

    for items in session.lists():
        for item in items:
            item.object = bpy.data.objects[item.name]

    bpy.context.scene.MNSession.load(session)

    # return database
