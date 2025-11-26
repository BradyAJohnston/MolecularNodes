import bpy
import databpy as db
import numpy as np
from bpy.types import Context, Depsgraph, Object
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


class GeometrySet:
    def __init__(self, obj: Object, context: None | Context = None):
        self.obj = obj
        self.context = context if isinstance(context, Context) else bpy.context
        depsgraph: Depsgraph = self.context.view_layer.depsgraph
        if depsgraph is None:
            raise ValueError

        self.eval_obj = depsgraph.id_eval_get(self.obj)
        self.geom = self.eval_obj.evaluated_geometry()

    @property
    def instances(self):
        return self.geom.instances_pointcloud()

    def instance_info(self):
        return {
            x: db.Attribute(self.instances.attributes[x])
            for x in ["instance_transform", ".reference_index"]
        }

    def mesh_info(self):
        mesh = self.geom.mesh
        # print("\n".join(dir(mesh)))
        string = "{} vertices, {} edges, {} polygons".format(
            len(mesh.vertices), len(mesh.edges), len(mesh.polygons)
        )
        return string

    def __repr__(self):
        string = "\n".join(
            [
                str(x)
                for x in [
                    "Mesh name: {}".format(self.geom.mesh.name),
                    "Mesh: {}".format(self.mesh_info()),
                    "Mesh Base: {}".format(self.geom.mesh_base),
                    "Pointcloud: {}".format(self.geom.pointcloud),
                    "Grease Pencil: {}".format(self.geom.grease_pencil),
                    "Curves: {}".format(self.geom.curves),
                    "Volumes: {}".format(self.geom.volume),
                    "Instances: {} points with {} unique values".format(
                        len(self.instance_info()["instance_transform"]),
                        len(np.unique(self.instance_info()[".reference_index"])),
                    ),
                ]
            ]
        )
        return string
