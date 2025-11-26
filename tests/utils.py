import bpy
import databpy as db
import numpy as np
from bpy.types import Context, Depsgraph, Object
from syrupy.extensions.amber import AmberSnapshotExtension


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

    def _format_attribute(self, attr: db.Attribute, max_display: int = 5) -> str:
        arr = attr.as_array()
        unique = np.unique(arr)
        n_unique = len(unique)

        parts = [f"shape={arr.shape}, atype={attr.atype.value}"]

        if n_unique == 1:
            parts.append(f"constant={unique[0]:.3f}")
        elif n_unique <= max_display:
            parts.append(f"unique={n_unique}, values={list(unique)}")
        else:
            parts.append(f"unique={n_unique}")
            if arr.dtype.kind in ["i", "f"]:
                parts.append(f"range=[{arr.min():.3g}, {arr.max():.3g}]")

        return ", ".join(parts)

    def _summarize_attributes(
        self, attributes_dict, label: str, max_attrs: int = 50
    ) -> list[str]:
        if not attributes_dict:
            return []

        lines = [f"\n{label}:"]
        attr_names = [
            name for name in attributes_dict.keys() if not name.startswith(".")
        ]
        attr_names.sort()

        for name in attr_names[:max_attrs]:
            attr = db.Attribute(attributes_dict[name])
            lines.append(f"  {name}: {self._format_attribute(attr)}")

        if len(attr_names) > max_attrs:
            lines.append(f"  ... and {len(attr_names) - max_attrs} more attributes")

        return lines

    def __repr__(self):
        lines = []

        mesh = self.geom.mesh
        if mesh:
            lines.append(f"Mesh: {mesh.name}")
            lines.append(
                f"  Geometry: {len(mesh.vertices)} verts, {len(mesh.edges)} edges, {len(mesh.polygons)} polys"
            )
            lines.extend(self._summarize_attributes(mesh.attributes, "  Attributes"))

        pointcloud = self.geom.pointcloud
        if pointcloud:
            if "position" in pointcloud.attributes:
                n_points = len(db.Attribute(pointcloud.attributes["position"]))
            else:
                n_points = 0
            lines.append(f"\nPointcloud: {n_points} points")
            lines.extend(
                self._summarize_attributes(pointcloud.attributes, "  Attributes")
            )

        instances = self.instances
        if instances:
            if "position" in instances.attributes:
                n_points = len(db.Attribute(instances.attributes["position"]))
            else:
                n_points = 0

            if n_points > 0:
                lines.append(f"\nInstances: {n_points} points")

                if "instance_transform" in instances.attributes:
                    transform_attr = db.Attribute(
                        instances.attributes["instance_transform"]
                    )
                    transforms = transform_attr.as_array()
                    lines.append(
                        f"  Transforms: shape={transforms.shape}, dtype={transforms.dtype.name}"
                    )

                    n_show = min(50, transforms.shape[0])
                    for i in range(n_show):
                        transform = transforms[i]
                        lines.append(f"    Transform[{i}]:")
                        for row_idx in range(4):
                            row = transform[row_idx]
                            row_str = f"      [{row[0]:7.3f}, {row[1]:7.3f}, {row[2]:7.3f}, {row[3]:7.3f}]"
                            lines.append(row_str)

                    if transforms.shape[0] > n_show:
                        lines.append(
                            f"    ... and {transforms.shape[0] - n_show} more transforms"
                        )

            if ".reference_index" in instances.attributes:
                ref_attr = db.Attribute(instances.attributes[".reference_index"])
                ref_arr = ref_attr.as_array()
                n_unique = len(np.unique(ref_arr))
                lines.append(f"  Unique instances: {n_unique}")
            lines.extend(
                self._summarize_attributes(instances.attributes, "  Attributes", 8)
            )

        if self.geom.curves:
            lines.append("\nCurves: present")

        if self.geom.volume:
            lines.append("\nVolume: present")

        if self.geom.grease_pencil:
            lines.append("\nGrease Pencil: present")

        return "\n".join(lines)
