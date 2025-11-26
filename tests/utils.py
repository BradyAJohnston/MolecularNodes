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

    def _get_point_count(self, attributes) -> int:
        if "position" in attributes:
            return len(db.Attribute(attributes["position"]))
        return 0

    def _format_value(self, value, dtype_kind: str) -> str:
        if dtype_kind == "b":
            return str(bool(value))
        elif dtype_kind == "i":
            return str(int(value))
        elif dtype_kind == "f":
            return f"{value:.3g}"
        else:
            return str(value)

    def _format_attribute(self, attr: db.Attribute, max_display: int = 5) -> str:
        arr = attr.as_array()
        dtype_kind = arr.dtype.kind

        if dtype_kind == "f":
            arr = np.where(np.isnan(arr) | np.isinf(arr), 0, arr)
            arr = np.where(np.abs(arr) < 1e-10, 0, arr)
            arr = np.where(np.abs(arr) > 1e10, 0, arr)

        unique = np.unique(arr)
        n_unique = len(unique)

        parts = [f"shape={arr.shape}, atype={attr.atype.value}"]

        if n_unique == 1:
            formatted_val = self._format_value(unique[0], dtype_kind)
            parts.append(f"constant={formatted_val}")
        elif n_unique <= max_display:
            formatted_values = [self._format_value(v, dtype_kind) for v in unique]
            parts.append(f"unique={n_unique}, values={formatted_values}")
        else:
            parts.append(f"unique={n_unique}")
            if dtype_kind in ["i", "f"]:
                min_val = self._format_value(arr.min(), dtype_kind)
                max_val = self._format_value(arr.max(), dtype_kind)
                parts.append(f"range=[{min_val}, {max_val}]")

        return ", ".join(parts)

    def _format_transform_matrix(
        self, transforms: np.ndarray, max_show: int = 50
    ) -> list[str]:
        lines = [
            f"  Transforms: shape={transforms.shape}, dtype={transforms.dtype.name}"
        ]

        n_show = min(max_show, transforms.shape[0])
        for i in range(n_show):
            transform = transforms[i]
            lines.append(f"    Transform[{i}]:")
            for row_idx in range(4):
                row = transform[row_idx]
                row_str = f"      [{row[0]:7.3f}, {row[1]:7.3f}, {row[2]:7.3f}, {row[3]:7.3f}]"
                lines.append(row_str)

        if transforms.shape[0] > n_show:
            lines.append(f"    ... and {transforms.shape[0] - n_show} more transforms")

        return lines

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

    def _summarize_mesh(self) -> list[str]:
        mesh = self.geom.mesh
        if not mesh:
            return []

        lines = [
            f"Mesh: {mesh.name}",
            f"  Geometry: {len(mesh.vertices)} verts, {len(mesh.edges)} edges, {len(mesh.polygons)} polys",
        ]
        lines.extend(self._summarize_attributes(mesh.attributes, "  Attributes"))
        return lines

    def _summarize_pointcloud(self) -> list[str]:
        pointcloud = self.geom.pointcloud
        if not pointcloud:
            return []

        n_points = self._get_point_count(pointcloud.attributes)
        lines = [f"\nPointcloud: {n_points} points"]
        lines.extend(self._summarize_attributes(pointcloud.attributes, "  Attributes"))
        return lines

    def _summarize_instances(self) -> list[str]:
        instances = self.instances
        if not instances:
            return []

        n_points = self._get_point_count(instances.attributes)
        if n_points == 0:
            return []

        lines = [f"\nInstances: {n_points} points"]

        if "instance_transform" in instances.attributes:
            transform_attr = db.Attribute(instances.attributes["instance_transform"])
            transforms = transform_attr.as_array()
            lines.extend(self._format_transform_matrix(transforms))

        if ".reference_index" in instances.attributes:
            ref_attr = db.Attribute(instances.attributes[".reference_index"])
            ref_arr = ref_attr.as_array()
            n_unique = len(np.unique(ref_arr))
            lines.append(f"  Unique instances: {n_unique}; {ref_arr}")
            lines.append(f"  Instance references: {self.geom.instance_references()}")

        lines.extend(
            self._summarize_attributes(instances.attributes, "  Attributes", 50)
        )
        return lines

    def _summarize_curves(self) -> list[str]:
        if self.geom.curves:
            return ["\nCurves: present"]
        return []

    def _summarize_volume(self) -> list[str]:
        if self.geom.volume:
            return ["\nVolume: present"]
        return []

    def _summarize_grease_pencil(self) -> list[str]:
        if self.geom.grease_pencil:
            return ["\nGrease Pencil: present"]
        return []

    def __repr__(self):
        lines = []
        lines.extend(self._summarize_mesh())
        lines.extend(self._summarize_pointcloud())
        lines.extend(self._summarize_instances())
        lines.extend(self._summarize_curves())
        lines.extend(self._summarize_volume())
        lines.extend(self._summarize_grease_pencil())
        return "\n".join(lines)
