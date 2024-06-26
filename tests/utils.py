import bpy

# from .conftest import molecularnodes as mn
import molecularnodes as mn

import numpy as np
import random
# import pathlib

from syrupy.extensions.amber import AmberSnapshotExtension


# we create a custom snapshot comparison class, which can handle numpy arrays
# and compare them properly. The class will serialize the numpy arrays into lists
# and when comparing them, reads the list back into a numpy array for comparison
# it checks for 'isclose' for floats and otherwise looks for absolute comparison
class NumpySnapshotExtension(AmberSnapshotExtension):
    def serialize(self, data, **kwargs):
        if isinstance(data, np.ndarray):
            return np.array2string(
                data, precision=1, threshold=1e3, floatmode="maxprec_equal"
            )
        return super().serialize(data, **kwargs)

    # def matches(self, *, serialized_data, snapshot_data):
    #     print(f"HELLOOOO")
    #     print(f"{serialized_data=}")
    #     print(f"{snapshot_data=}")
    #     serialized_data = np.array(ast.literal_eval(serialized_data)),
    #     snapshot_data = np.array(ast.literal_eval(snapshot_data)),
    #     print(f"{serialized_data=}")
    #     print(f"{snapshot_data=}")

    #     # super().assert_match(snapshot_custom, test_value)
    #     # def assert_match(self, snapshot_custom, test_value):
    #     if isinstance(serialized_data, np.ndarray):
    #         # if the values are floats, then we use a rough "isclose" to compare them
    #         # which helps with floating point issues. Between platforms geometry nodes
    #         # outputs some differences in the meshes which are usually off by ~0.01 or so

    #         else:
    #             assert (serialized_data == np.array(snapshot_data)).all()

    # else:
    #     super().matches(serialized_data=serialized_data, snapshot_data=snapshot_data)


def sample_attribute(
    object, attribute, n=100, evaluate=True, error: bool = False, seed=6
):
    if isinstance(object, mn.io.molecule.Molecule):
        object = object.object

    random.seed(seed)
    if error:
        attribute = mn.blender.obj.get_attribute(object, attribute, evaluate=evaluate)
        length = len(attribute)

        if n > length:
            idx = range(length)
        else:
            idx = random.sample(range(length), n)

        if len(attribute.data.shape) == 1:
            return attribute[idx]

        return attribute[idx, :]
    else:
        try:
            attribute = mn.blender.obj.get_attribute(
                object=object, name=attribute, evaluate=evaluate
            )
            length = len(attribute)

            if n > length:
                idx = range(length)
            else:
                idx = random.sample(range(length), n)

            if len(attribute.data.shape) == 1:
                return attribute[idx]

            return attribute[idx, :]
        except AttributeError as e:
            return np.array(e)


def sample_attribute_to_string(
    object, attribute, n=100, evaluate=True, precision=3, seed=6
):
    if isinstance(object, mn.io.molecule.Molecule):
        object = object.object
    try:
        array = sample_attribute(
            object, attribute=attribute, n=n, evaluate=evaluate, seed=seed
        )
    except AttributeError as e:
        print(f"Error {e}, unable to sample attribute {attribute} from {object}")
        return str(e)

    if array.dtype != bool:
        array = np.round(array, precision)
    length = len(array)
    threshold = 4 * length

    if n > length:
        idx = range(length)
    else:
        idx = random.sample(range(length), n)

    dimensions = len(np.shape(attribute))

    if dimensions == 1:
        array = attribute[idx]
    elif dimensions == 2:
        array = attribute[idx, :]

    return np.array2string(array, precision=precision, threshold=threshold)


def remove_all_molecule_objects(mda_session):
    for object in bpy.data.objects:
        try:
            obj_type = object["type"]
            if obj_type == "molecule":
                bpy.data.objects.remove(object)
        except KeyError:
            pass
    # remove frame change
    bpy.context.scene.frame_set(0)

    mda_session.universe_reps = {}
    mda_session.atom_reps = {}
    mda_session.rep_names = []
