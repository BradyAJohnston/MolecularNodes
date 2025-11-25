import bpy
import numpy as np
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


# https://docs.blender.org/api/4.5/bpy.types.GeometrySet.html


def get_geometry_set(
    obj: bpy.types.Object, context: None | bpy.types.Context = None
) -> bpy.types.GeometrySet:  # type: ignore
    if not context:
        context = bpy.context

    depsgraph = context.view_layer.depsgraph  # type: ignore
    if depsgraph is None:
        raise ValueError

    ob_eval = depsgraph.id_eval_get(obj)
    geom: bpy.types.GeometrySet = ob_eval.evaluated_geometry()  # type: ignore

    return geom


def geometry_set_to_dict(geom_set: bpy.types.GeometrySet) -> dict:
    return {}
