import pickle
import bpy
from bpy.app.handlers import persistent
from pathlib import Path


def stashpath():
    return f"{bpy.data.filepath}.MN_database"


@persistent
def _stash_save(scene: bpy.types.Scene) -> None:
    path = stashpath()
    stash_database(database=bpy.context.scene.MN_database, filepath=stashpath())
    print(f"Save database to: {path}")


@persistent
def _stash_load(scene: bpy.types.Context) -> None:
    apply_database(stashpath())


def stash_database(database, filepath) -> None:
    for item in database:
        if isinstance(item.object, bpy.types.Object):
            item.object = item.object.name
    with open(stashpath(), "wb") as f:
        pickle.dump(database, f)


def apply_database(file: Path) -> None:
    with open(file, "rb") as f:
        database = pickle.load(f)

    for item in database:
        item.object = bpy.data.objects[item.object]

        bpy.context.scene.MN_database.append(item)

    # return database
