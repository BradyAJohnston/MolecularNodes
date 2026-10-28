import bpy
from bpy.types import Collection
from databpy.collection import create_collection


def mn() -> Collection:
    "Return the 'MolecularNodes' collection, creating it first if required"
    return create_collection("Molecular Nodes")


def data() -> Collection:
    "Return the MolecularNodes/data collection and disable it"
    name = ".MN_data"

    try:
        return bpy.data.collections[name]
    except KeyError:
        collection = create_collection(name=name, parent=mn())

        assert bpy.context.view_layer is not None
        bpy.context.view_layer.layer_collection.children["Molecular Nodes"].children[
            collection.name
        ].exclude = True
        return collection


def frames(name: str = "") -> Collection:
    "Return a collection for storing the objects that are the frames of a trajectory"
    return create_collection(f".data_{name}_frames", parent=data())


def cellpack(name: str = "") -> Collection:
    "Return a new collection for storing the instances for a CellPack Ensemble"
    # always a new collection: reusing one by name would mix the instances of
    # separate imports
    collection = bpy.data.collections.new(f"cellpack_{name}")
    data().children.link(collection)
    return collection
