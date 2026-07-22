import databpy as db
import numpy as np
from syrupy.extensions.amber import AmberSnapshotExtension

# Number of significant digits that float values are rounded to before any summary
# statistic is derived from them. Blender's geometry nodes evaluate on the CPU with
# platform-specific maths, so the last bits of a float are not portable - snapshots
# have to compare quantities that survive that.
FLOAT_SIG_DIGITS = 3

# Counts of procedurally generated geometry (surfaces, ribbons, instanced meshes) drift
# by ~0.01% between the Linux/macOS/Windows Blender builds, so they are reported to this
# many significant digits rather than exactly.
COUNT_SIG_DIGITS = 3


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


def round_significant(arr: np.ndarray, digits: int = FLOAT_SIG_DIGITS) -> np.ndarray:
    """
    Round a float array to a number of significant digits.

    Rounding relative to each value's magnitude rather than to a fixed number of decimal
    places keeps small-magnitude attributes meaningful, while collapsing the last-bit
    differences that different platforms produce for the same computation.
    """
    arr = np.asarray(arr, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        magnitude = np.where(arr == 0, 0, np.floor(np.log10(np.abs(arr))))
    factor = 10.0 ** (digits - 1 - magnitude)
    return np.where(arr == 0, 0.0, np.round(arr * factor) / factor)


class GeometrySet(db.GeometrySet):
    """
    `databpy.GeometrySet`, extended with attribute summaries for snapshot testing.

    Call `summary()` to get a dict for snapshotting. Use `strict=True` for geometry that
    comes straight from a parsed file, where every count is bit-exact on every platform,
    and the default `strict=False` for geometry that has been through geometry nodes,
    where counts are quantised so that platform-level meshing differences don't fail the
    snapshot.
    """

    def __init__(self, obj, context=None, strict: bool = False):
        super().__init__(obj, context)
        self.strict = strict

    def _count(self, n: int) -> str:
        """Format a count, quantising generated geometry so small drift doesn't show."""
        if self.strict or n < 10**COUNT_SIG_DIGITS:
            return str(n)
        step = 10 ** (len(str(n)) - COUNT_SIG_DIGITS)
        return f"~{round(n / step) * step}"

    def _shape(self, shape: tuple[int, ...]) -> str:
        dims = ", ".join(self._count(dim) for dim in shape)
        return f"({dims},)" if len(shape) == 1 else f"({dims})"

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
            return f"{value:.{FLOAT_SIG_DIGITS}g}"
        else:
            return str(value)

    def _format_attribute(self, attr: db.Attribute, max_display: int = 5) -> str:
        arr = attr.as_array()
        dtype_kind = arr.dtype.kind

        if dtype_kind == "f":
            arr = np.where(np.isnan(arr) | np.isinf(arr), 0, arr)
            arr = np.where(np.abs(arr) < 1e-8, 0, arr)
            arr = np.where(np.abs(arr) > 1e10, 0, arr)
            # quantise before deriving anything, so that unique counts and min/max are
            # not sensitive to the last bits of the float
            arr = round_significant(arr)

        unique = np.unique(arr)
        n_unique = len(unique)

        parts = [f"shape={self._shape(arr.shape)}, atype={attr.atype.value}"]

        if n_unique == 1:
            parts.append(f"constant={self._format_value(unique[0], dtype_kind)}")
        elif n_unique <= max_display:
            values = [self._format_value(v, dtype_kind) for v in unique]
            parts.append(f"unique={n_unique}, values={values}")
        else:
            # a bare cardinality is only reported for exactly-reproducible geometry -
            # a couple of extra vertices from a mesher changes it without meaning
            # anything, and floats have no stable cardinality at all
            if self.strict and dtype_kind in ["i", "b"]:
                parts.append(f"unique={n_unique}")
            if dtype_kind in ["i", "f"]:
                min_val = self._format_value(arr.min(), dtype_kind)
                max_val = self._format_value(arr.max(), dtype_kind)
                parts.append(f"range=[{min_val}, {max_val}]")

        return ", ".join(parts)

    def _format_transforms(self, transforms: np.ndarray, max_show: int = 50) -> list:
        rounded = np.round(transforms, 3)
        # -0.0 and 0.0 format differently but are the same value, and which one a
        # computation lands on is not portable
        rounded[rounded == 0] = 0.0
        lines = [
            " ".join(
                "[" + ", ".join(f"{value:.3f}" for value in row) + "]"
                for row in rounded[i]
            )
            for i in range(min(max_show, transforms.shape[0]))
        ]
        if transforms.shape[0] > max_show:
            lines.append(f"... and {transforms.shape[0] - max_show} more transforms")
        return lines

    def _format_reference(self, reference) -> str:
        """
        Describe an instance reference without its exact topology.

        The default repr of a referenced geometry embeds vertex/edge/face counts, which
        for generated geometry are not portable between platforms.
        """
        name = getattr(reference, "name", None)
        if name:
            return f"{type(reference).__name__}({name})"

        mesh = getattr(reference, "mesh", None)
        if mesh is not None:
            return (
                f"{type(reference).__name__}("
                f"{self._count(len(mesh.vertices))} verts, "
                f"{self._count(len(mesh.polygons))} polys)"
            )
        return type(reference).__name__

    def _summarize_attributes(self, attributes_dict, max_attrs: int = 50) -> dict:
        if not attributes_dict:
            return {}

        attr_names = sorted(
            name for name in attributes_dict.keys() if not name.startswith(".")
        )

        summary = {
            name: self._format_attribute(db.Attribute(attributes_dict[name]))
            for name in attr_names[:max_attrs]
        }
        if len(attr_names) > max_attrs:
            summary["..."] = f"and {len(attr_names) - max_attrs} more attributes"

        return summary

    def _summarize_mesh(self) -> dict:
        mesh = self.mesh
        if not mesh:
            return {}

        return {
            "mesh": {
                # `mesh.name` is the evaluated depsgraph mesh's data-block name, which
                # Blender does not guarantee to be stable (it can be a generic "Mesh"
                # depending on version/build) - use the object's own name instead.
                "name": self.object.name,
                "verts": self._count(len(mesh.vertices)),
                "edges": self._count(len(mesh.edges)),
                "polys": self._count(len(mesh.polygons)),
                "attributes": self._summarize_attributes(mesh.attributes),
            }
        }

    def _summarize_pointcloud(self) -> dict:
        pointcloud = self.pointcloud
        if not pointcloud:
            return {}

        return {
            "pointcloud": {
                "points": self._count(self._get_point_count(pointcloud.attributes)),
                "attributes": self._summarize_attributes(pointcloud.attributes),
            }
        }

    def _summarize_instances(self) -> dict:
        instances = self.instances
        if not instances:
            return {}

        n_points = self._get_point_count(instances.attributes)
        if n_points == 0:
            return {}

        summary = {"points": self._count(n_points)}

        if "instance_transform" in instances.attributes:
            transform_attr = db.Attribute(instances.attributes["instance_transform"])
            summary["transforms"] = self._format_transforms(transform_attr.as_array())

        if ".reference_index" in instances.attributes:
            ref_arr = db.Attribute(instances.attributes[".reference_index"]).as_array()
            summary["unique_references"] = len(np.unique(ref_arr))
            summary["reference_index"] = str(ref_arr)
            summary["references"] = [
                self._format_reference(reference)
                for reference in self.instance_references
            ]

        filtered_attrs = {
            k: v for k, v in instances.attributes.items() if k != "position"
        }
        summary["attributes"] = self._summarize_attributes(filtered_attrs)
        return {"instances": summary}

    def summary(self) -> dict:
        """A snapshot-friendly summary of every component of the evaluated geometry."""
        summary: dict = {}
        summary.update(self._summarize_mesh())
        summary.update(self._summarize_pointcloud())
        summary.update(self._summarize_instances())
        if self.curves:
            summary["curves"] = "present"
        if self.volume:
            summary["volume"] = "present"
        if self.grease_pencil:
            summary["grease_pencil"] = "present"
        return summary
