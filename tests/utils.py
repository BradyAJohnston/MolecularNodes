import io
import os
from pathlib import Path
from typing import Any
import bpy
import databpy as db
import numpy as np
from syrupy.extensions.amber import AmberSnapshotExtension
from syrupy.extensions.single_file import SingleFileSnapshotExtension

# Number of significant digits that float values are rounded to before any summary
# statistic is derived from them. Blender's geometry nodes evaluate on the CPU with
# platform-specific maths, so the last bits of a float are not portable - snapshots
# have to compare quantities that survive that.
FLOAT_SIG_DIGITS = 3

# Counts of procedurally generated geometry (surfaces, ribbons, instanced meshes) drift
# by ~0.01% between the Linux/macOS/Windows Blender builds, so they are reported to this
# many significant digits rather than exactly.
COUNT_SIG_DIGITS = 3


#: the operators Molecular Nodes registers under `bpy.ops.mn`, which the bpy stubs
#: can't know about
mn_ops: Any = bpy.ops.mn  # ty: ignore[unresolved-attribute]


def scene() -> bpy.types.Scene:
    """The current scene, which the bpy stubs type as optional."""
    assert bpy.context.scene is not None
    return bpy.context.scene


def view_layer() -> bpy.types.ViewLayer:
    """The current view layer, which the bpy stubs type as optional."""
    assert bpy.context.view_layer is not None
    return bpy.context.view_layer


def node_input(node: bpy.types.Node, key: str | int) -> bpy.types.NodeSocket:
    """An input socket of a node, whose inputs the stubs type as optional."""
    assert node.inputs is not None
    return node.inputs[key]


def socket_links(socket: bpy.types.NodeSocket) -> "list[bpy.types.NodeLink]":
    """The links of a node socket, which the stubs type as optional."""
    assert socket.links is not None
    return list(socket.links)


def object_mesh(obj: bpy.types.Object) -> bpy.types.Mesh:
    """The data of an object that must be a mesh."""
    assert isinstance(obj.data, bpy.types.Mesh)
    return obj.data


def interface_items(
    tree: bpy.types.NodeTree,
) -> "bpy.types.bpy_prop_collection[bpy.types.NodeTreeInterfaceItem]":
    """The interface items of a node tree, whose interface the stubs type as optional."""
    assert tree.interface is not None
    return tree.interface.items_tree


def interface_socket(
    tree: bpy.types.NodeTree, key: str | int
) -> bpy.types.NodeTreeInterfaceSocket:
    """An interface item of a node tree that must be a socket."""
    item = interface_items(tree)[key]
    assert isinstance(item, bpy.types.NodeTreeInterfaceSocket)
    return item


def modifier_node_group(obj: bpy.types.Object) -> bpy.types.NodeTree:
    """The node group of an object's first modifier, which must be a nodes modifier."""
    modifier = obj.modifiers[0]
    assert isinstance(modifier, bpy.types.NodesModifier)
    assert modifier.node_group is not None
    return modifier.node_group


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
                data, precision=1, threshold=2000, floatmode="maxprec_equal"
            )
        return super().serialize(data, **kwargs)


# Tolerances for golden-image comparison, following Blender's own render tests
# (tests/python/modules/render_report.py): a pixel fails when any channel
# differs by more than FAIL_THRESHOLD, and the image fails when more than
# FAIL_PERCENT of pixels fail. Blender uses 0.016 (~4/255) and 1% for Cycles.
IMAGE_FAIL_THRESHOLD = 4 / 255
IMAGE_FAIL_PERCENT = 1.0

# where a failing comparison writes the received image and the difference
# image for inspection (uploadable as a CI artifact)
IMAGE_FAILURE_DIR = Path(__file__).parent / "image_failures"


def _decode_png(data: bytes) -> np.ndarray:
    from PIL import Image as PILImage

    return np.asarray(PILImage.open(io.BytesIO(data)).convert("RGBA"), dtype=np.uint8)


def compare_images(
    reference: bytes,
    received: bytes,
    fail_threshold: float = IMAGE_FAIL_THRESHOLD,
    fail_percent: float = IMAGE_FAIL_PERCENT,
) -> tuple[bool, str]:
    """
    Compare two PNGs with Blender-style tolerances.

    Returns (matches, stats) where stats describes the max channel error and
    the percentage of failing pixels.
    """
    ref = _decode_png(reference)
    new = _decode_png(received)
    if ref.shape != new.shape:
        return False, f"shape mismatch: {ref.shape} != {new.shape}"
    diff = np.abs(ref.astype(np.float32) - new.astype(np.float32)) / 255
    failing = (diff > fail_threshold).any(axis=-1)
    percent = failing.mean() * 100
    stats = f"max error = {diff.max():.4f}, {percent:.2f}% of pixels over {fail_threshold:.4f} (allowed {fail_percent}%)"
    return bool(percent <= fail_percent), stats


class ImageSnapshotExtension(SingleFileSnapshotExtension):
    """
    Golden-image snapshots: `.png` files compared with a tolerance.

    Bit-exact renders are not portable across platforms, so snapshots match
    when the difference stays within Blender's own render-test tolerances.
    Update the goldens with `pytest --snapshot-update`; on failure the
    received image and a difference image are written to
    ``tests/image_failures/`` for inspection.
    """

    file_extension = "png"
    fail_threshold = IMAGE_FAIL_THRESHOLD
    fail_percent = IMAGE_FAIL_PERCENT

    @classmethod
    def with_tolerance(
        cls,
        fail_threshold: float = IMAGE_FAIL_THRESHOLD,
        fail_percent: float = IMAGE_FAIL_PERCENT,
    ) -> type["ImageSnapshotExtension"]:
        """A variant of this extension with looser comparison tolerances."""
        return type(
            cls.__name__,
            (cls,),
            {"fail_threshold": fail_threshold, "fail_percent": fail_percent},
        )

    def matches(self, *, serialized_data, snapshot_data) -> bool:
        try:
            ok, stats = compare_images(
                bytes(snapshot_data),
                bytes(serialized_data),
                fail_threshold=self.fail_threshold,
                fail_percent=self.fail_percent,
            )
        except Exception as e:
            return self._dump_failure(serialized_data, f"decode failed: {e}")
        if not ok:
            return self._dump_failure(serialized_data, stats, bytes(snapshot_data))
        return True

    def _dump_failure(
        self, received, stats: str, reference: bytes | None = None
    ) -> bool:
        self._failure_stats = stats
        IMAGE_FAILURE_DIR.mkdir(exist_ok=True)
        # matches() has no test context, but pytest exposes the running test
        current_test = os.environ.get("PYTEST_CURRENT_TEST", "image (call)")
        name = current_test.split("::")[-1].split(" ")[0]
        received_file = IMAGE_FAILURE_DIR / f"{name}.received.png"
        received_file.write_bytes(bytes(received))
        self._failure_files = [received_file]
        if reference is not None:
            from PIL import Image as PILImage

            diff = np.abs(
                _decode_png(reference).astype(np.int16)
                - _decode_png(bytes(received)).astype(np.int16)
            )
            # amplify so small differences are visible, and make it opaque
            diff = np.clip(diff * 4, 0, 255).astype(np.uint8)
            diff[..., 3] = 255
            diff_file = IMAGE_FAILURE_DIR / f"{name}.diff.png"
            PILImage.fromarray(diff).save(diff_file)
            self._failure_files.append(diff_file)
        return False

    def diff_lines(self, serialized_data, snapshot_data):
        # a byte-level diff of two PNGs is noise - report the stats and where
        # the images for inspection were written instead
        yield getattr(self, "_failure_stats", "images differ")
        for file in getattr(self, "_failure_files", []):
            yield f"see {file}"


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

        if dtype_kind == "b":
            # counts rather than the unique values, so a change in how many points
            # are flagged shows up in the snapshot
            n_true = int(np.count_nonzero(arr))
            parts.append(
                f"true={self._count(n_true)}, false={self._count(arr.size - n_true)}"
            )
        elif n_unique == 1:
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


def sphere_value(mol) -> str | None:
    "The `Sphere` value on whichever style node has one."
    for node in mol.modifier_node_tree.nodes:
        for socket in node.inputs:
            if socket.name == "Sphere":
                return socket.default_value
    return None
