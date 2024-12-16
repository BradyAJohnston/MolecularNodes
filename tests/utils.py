import bpy
import molecularnodes as mn

import numpy as np
import random

from syrupy.extensions.amber import AmberSnapshotExtension


# we create a custom snapshot comparison class, which can handle numpy arrays
# and compare them properly. The class will serialize the numpy arrays into lists
# and when comparing them, reads the list back into a numpy array for comparison
# it checks for 'isclose' for floats and otherwise looks for absolute comparison
class NumpySnapshotExtension(AmberSnapshotExtension):
    def __init__(self):
        super().__init__()
        self.custom_suffix: str | None = None

    def serialize(self, data, cutoff=1000, **kwargs):
        
        if isinstance(data, np.ndarray):
            shape = data.shape
            if len(shape) == 1:
                if len(data) > cutoff:
                    data = data[:cutoff]
                else:
                    data = data[: int(cutoff / 10),]

            return np.array2string(
                data, precision=1, threshold=2e3, floatmode="maxprec_equal"
            )
        return super().serialize(data, **kwargs)

    

def sample_attribute(
    object, attribute, n=100, evaluate=True, error: bool = False, seed=6
):
    if isinstance(object, mn.entities.molecule.Molecule):
        object = object.object

    random.seed(seed)
    if error:
        attribute = mn.bpyd.named_attribute(
            obj=object, name=attribute, evaluate=evaluate
        )
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
            attribute = mn.bpyd.named_attribute(
                obj=object, name=attribute, evaluate=evaluate
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
    if isinstance(object, mn.entities.molecule.Molecule):
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
