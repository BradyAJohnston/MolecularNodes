"""
Import a MolViewSpec ``.mvsj`` state as Molecular Nodes entities.

The MVS tree is walked once, inheriting context (download url, parse format,
structure params) down each branch. Each ``structure`` node becomes one
:class:`Molecule`; each ``representation`` under a ``component`` becomes one
``add_style`` branch, with selections and colors resolved to named attributes
(see ``docs/dev/mvs-import.md`` for the overall design).

Anything the importer cannot faithfully translate is collected and reported
as a single warning at the end - a valid file always yields a best-effort
scene rather than an error.
"""

import hashlib
import urllib.parse
import urllib.request
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import bpy
import numpy as np
from databpy import AttributeTypes
from matplotlib.colors import to_rgba
from ...download import CACHE_DIR
from ..molecule.base import Molecule
from .selectors import component_mask

# MVS representation type -> MN style name; the second group are nearest-style
# approximations that import with a warning
STYLE_MAPPING = {
    "cartoon": "cartoon",
    "ball_and_stick": "ball_and_stick",
    "spacefill": "spheres",
    "surface": "surface",
}
STYLE_FALLBACKS = {
    "putty": "ribbon",
    "backbone": "ribbon",
    "line": "sticks",
    "carbohydrate": "ball_and_stick",
}

PARSE_SUFFIXES = {"mmcif": ".cif", "bcif": ".bcif", "pdb": ".pdb"}

# MVS coordinates and distances are in Angstrom, the Blender world in nm
_WORLD_SCALE = 0.1

# node kinds that are recognised but out of scope for the current phase
_UNSUPPORTED_KINDS = {
    "label",
    "label_from_source",
    "label_from_uri",
    "tooltip",
    "tooltip_from_source",
    "tooltip_from_uri",
    "primitives",
    "primitives_from_uri",
    "primitive",
    "transform",
    "instance",
    "volume",
    "volume_representation",
    "clip",
    "mesh_from_source",
    "mesh_from_uri",
    "component_from_source",
    "component_from_uri",
    "color_from_source",
    "color_from_uri",
    "coordinates",
    "shape",
    "transition",
}


# kinds handled in their expected context, reported differently when reached
# under an unsupported parent node
_CONTEXT_KINDS = {
    "color",
    "opacity",
    "representation",
    "component",
    "focus",
    "structure",
    "download",
    "parse",
}


class MVSImportWarning(UserWarning):
    """Parts of a MolViewSpec file could not be faithfully imported."""


@dataclass
class _Context:
    """Inherited state while walking down the MVS tree."""

    url: str | None = None
    format: str | None = None


@dataclass
class _SceneState:
    """Scene-level results collected during the walk."""

    background: tuple[float, float, float, float] | None = None
    camera: dict[str, Any] | None = None
    # positions (world-space) and node params of the most recent focused component
    focus_points: np.ndarray | None = None
    focus_params: dict[str, Any] | None = None


def load(file_path: str | Path) -> list[Molecule]:
    """
    Import a MolViewSpec ``.mvsj`` file.

    Each ``structure`` in the state becomes a :class:`Molecule` with one style
    branch per representation. Scene-level ``canvas``, ``camera`` and
    ``focus`` nodes are applied to the current Blender scene. Unsupported
    parts of the file are skipped with a single summary warning.

    Parameters
    ----------
    file_path : str | Path
        Path to the ``.mvsj`` file. Download URLs inside the file may be
        absolute ``http(s)`` URLs (cached on disk) or paths relative to the
        file itself.

    Returns
    -------
    list[Molecule]
        The imported entities, in file order.
    """
    return MVSImporter(file_path).import_state()


class MVSImporter:
    def __init__(self, file_path: str | Path) -> None:
        self.file_path = Path(file_path)
        self.base_dir = self.file_path.parent
        self.warnings: list[str] = []
        self.molecules: list[Molecule] = []
        self.scene_state = _SceneState()
        self._counter = 0
        self.root = self._load_root()

    def _load_root(self):
        # imported lazily so the addon still loads if molviewspec is missing
        from molviewspec import MVSJ

        data = MVSJ.load(self.file_path).data
        if getattr(data, "kind", "single") == "multiple":
            snapshots = data.snapshots
            if len(snapshots) > 1:
                self.warnings.append(
                    f"multi-state files are not yet supported; importing the "
                    f"first of {len(snapshots)} snapshots"
                )
            return snapshots[0].root
        return data.root

    # ------------------------------------------------------------------ walk
    def import_state(self) -> list[Molecule]:
        self._walk(self.root, _Context())
        self._apply_scene_state()
        if self.warnings:
            unique = list(dict.fromkeys(self.warnings))
            warnings.warn(
                "Parts of the MolViewSpec file were not imported:\n- "
                + "\n- ".join(unique),
                category=MVSImportWarning,
                stacklevel=2,
            )
        return self.molecules

    def _walk(self, node, context: _Context) -> None:
        for child in node.children or []:
            params = dict(child.params or {})
            match child.kind:
                case "download":
                    self._walk(child, _Context(url=params.get("url")))
                case "parse":
                    self._walk(
                        child,
                        _Context(url=context.url, format=params.get("format")),
                    )
                case "structure":
                    self._import_structure(child, params, context)
                case "canvas":
                    background = params.get("background_color")
                    if background is not None:
                        self.scene_state.background = to_rgba(background)
                    self._walk(child, context)
                case "camera":
                    self.scene_state.camera = params
                case _:
                    self._warn_unsupported(child.kind)
                    self._walk(child, context)

    def _warn_unsupported(self, kind: str) -> None:
        if kind in _UNSUPPORTED_KINDS:
            self.warnings.append(f"'{kind}' nodes are not yet supported")
        elif kind in _CONTEXT_KINDS:
            # a supported kind reached under an unsupported parent, e.g. a
            # color node inside a volume representation
            self.warnings.append(
                f"'{kind}' nodes under unsupported parents are skipped"
            )
        else:
            self.warnings.append(f"unknown node kind '{kind}'")

    # ------------------------------------------------------------- structure
    def _import_structure(self, node, params: dict, context: _Context) -> None:
        if context.url is None or context.format is None:
            self.warnings.append(
                "structure node without a download/parse context; skipping it"
            )
            return
        if context.format not in PARSE_SUFFIXES:
            self.warnings.append(
                f"parse format '{context.format}' is not supported; skipping "
                "the structure"
            )
            return

        source = self._resolve_source(context.url, context.format)
        if source is None:
            return

        name = Path(urllib.parse.urlparse(context.url).path).stem or "MVS"
        mol = Molecule.load(source, name=name)
        self.molecules.append(mol)

        structure_type = params.get("type", "model")
        assembly = structure_type == "assembly"
        if structure_type in ("symmetry", "symmetry_mates"):
            self.warnings.append(
                f"structure type '{structure_type}' is not supported; importing "
                "the deposited model instead"
            )
        if assembly and params.get("assembly_id") not in (None, "1"):
            self.warnings.append(
                f"only the first assembly is supported; requested assembly "
                f"'{params.get('assembly_id')}'"
            )
        for key in ("model_index", "block_index", "block_header"):
            if params.get(key) not in (None, 0):
                self.warnings.append(f"structure parameter '{key}' is ignored")

        label_fields = self._label_fields(mol, source)

        for child in node.children or []:
            child_params = dict(child.params or {})
            match child.kind:
                case "component":
                    self._import_component(
                        child, child_params, mol, label_fields, assembly
                    )
                case _:
                    self._warn_unsupported(child.kind)

    def _resolve_source(self, url: str, format: str) -> Path | None:
        """A local file path for a download node's url, fetching if remote."""
        suffix = PARSE_SUFFIXES[format]
        parsed = urllib.parse.urlparse(url)
        if parsed.scheme in ("http", "https"):
            cache_dir = CACHE_DIR / "mvs"
            cache_dir.mkdir(parents=True, exist_ok=True)
            name = hashlib.sha256(url.encode()).hexdigest()[:16] + suffix
            target = cache_dir / name
            if not target.exists():
                try:
                    urllib.request.urlretrieve(url, target)
                except Exception as error:
                    self.warnings.append(f"failed to download '{url}': {error}")
                    return None
            return target
        if parsed.scheme == "file":
            return Path(urllib.request.url2pathname(parsed.path))
        path = Path(url)
        if not path.is_absolute():
            path = self.base_dir / path
        if not path.exists():
            self.warnings.append(f"source file '{url}' not found; skipping")
            return None
        return path

    def _label_fields(self, mol: Molecule, source: Path) -> dict[str, np.ndarray]:
        """
        ``label_asym_id`` / ``label_seq_id`` per atom, from a secondary parse.

        The main import keeps the author-assigned identifiers, but MVS
        expressions may select on the standard (label) ones. Only mmCIF-family
        files carry them; other formats fall back to auth values downstream.
        """
        if source.suffix not in (".cif", ".bcif"):
            return {}
        try:
            from biotite.structure.io import pdbx

            file_cls = pdbx.CIFFile if source.suffix == ".cif" else pdbx.BinaryCIFFile
            array = pdbx.get_structure(
                file_cls.read(source),
                model=1,
                extra_fields=["label_asym_id", "label_seq_id"],
            )
            if array.array_length() != mol.universe.atoms.n_atoms:
                return {}
            # hetero atoms carry "." instead of a number, which maps to -1 so
            # it can never match a requested sequence position
            raw_seq = array.label_seq_id.astype(str)
            seq_ids = np.full(len(raw_seq), -1, dtype=int)
            numeric = np.char.isdigit(raw_seq)
            seq_ids[numeric] = raw_seq[numeric].astype(int)
            return {
                "label_asym_id": array.label_asym_id.astype(str),
                "label_seq_id": seq_ids,
            }
        except Exception:
            return {}

    # ------------------------------------------------------------- component
    def _import_component(
        self,
        node,
        params: dict,
        mol: Molecule,
        label_fields: dict[str, np.ndarray],
        assembly: bool,
    ) -> None:
        mask = component_mask(mol, params.get("selector"), label_fields, self.warnings)

        selection_name = None
        if not mask.all():
            selection_name = self._store_attribute(mol, mask, "selection")

        representations = []
        for child in node.children or []:
            child_params = dict(child.params or {})
            match child.kind:
                case "representation":
                    representations.append((child, child_params))
                case "focus":
                    positions = mol.named_attribute("position")[mask]
                    if len(positions):
                        self.scene_state.focus_points = positions
                        self.scene_state.focus_params = child_params
                case _:
                    self._warn_unsupported(child.kind)

        for child, child_params in representations:
            self._import_representation(
                child, child_params, mol, mask, selection_name, label_fields, assembly
            )

    def _store_attribute(self, mol: Molecule, data: np.ndarray, tag: str) -> str:
        self._counter += 1
        name = f"mvs_{tag}_{self._counter}"
        atype = (
            AttributeTypes.BOOLEAN if data.dtype == bool else AttributeTypes.FLOAT_COLOR
        )
        mol.store_named_attribute(data, name=name, atype=atype)
        return name

    # -------------------------------------------------------- representation
    def _import_representation(
        self,
        node,
        params: dict,
        mol: Molecule,
        mask: np.ndarray,
        selection_name: str | None,
        label_fields: dict[str, np.ndarray],
        assembly: bool,
    ) -> None:
        mvs_type = params.get("type", "cartoon")
        style = STYLE_MAPPING.get(mvs_type)
        if style is None:
            style = STYLE_FALLBACKS.get(mvs_type)
            if style is None:
                self.warnings.append(
                    f"representation '{mvs_type}' is not supported; skipping it"
                )
                return
            self.warnings.append(
                f"representation '{mvs_type}' is approximated with the '{style}' style"
            )

        colors: list[tuple[str, Any]] = []
        opacity = 1.0
        for child in node.children or []:
            child_params = dict(child.params or {})
            match child.kind:
                case "color":
                    if child_params.get("color") is not None:
                        colors.append(
                            (child_params["color"], child_params.get("selector"))
                        )
                case "opacity":
                    opacity *= float(child_params.get("opacity", 1.0))
                case _:
                    self._warn_unsupported(child.kind)

        color_input = None
        if colors or opacity != 1.0:
            rgba = np.asarray(mol.named_attribute("Color"), dtype=np.float32).copy()
            for color_value, sub_selector in colors:
                sub_mask = mask.copy()
                if sub_selector is not None:
                    sub_mask &= component_mask(
                        mol, sub_selector, label_fields, self.warnings
                    )
                rgba[sub_mask, :3] = to_rgba(color_value)[:3]
            # the default materials read the shader alpha from the color
            # attribute's alpha channel
            rgba[mask, 3] = opacity
            color_input = self._store_attribute(mol, rgba, "color")

        mol.add_style(
            style=style,
            selection=selection_name,
            color=color_input,
            assembly=assembly,
        )

    # ------------------------------------------------------------- the scene
    def _apply_scene_state(self) -> None:
        state = self.scene_state
        if (
            state.background is None
            and state.camera is None
            and (state.focus_points is None)
        ):
            return

        # applied through the standalone scene helpers rather than a Canvas,
        # which would reset render settings the user (or a test) has chosen
        from ...scene.camera import Camera
        from ...scene.world import WorldTree

        if bpy.context.scene.camera is None:
            camera_data = bpy.data.cameras.new("Camera")
            camera_object = bpy.data.objects.new("Camera", camera_data)
            bpy.context.scene.collection.objects.link(camera_object)
            bpy.context.scene.camera = camera_object

        if state.background is not None:
            world = WorldTree(bpy.context.scene)
            try:
                world.background = state.background
            except ValueError:
                # a bare scene without the MN world shader: build a minimal one
                from nodebpy import shader

                with world.reset() as surface:
                    shader.Background(color=state.background) >> surface
        camera = Camera()
        if state.camera is not None:
            self._apply_camera(camera, state.camera)
        elif state.focus_points is not None:
            # MVS focus defines the view direction (default looking down -Z);
            # orient the camera first, then solve its position over the points
            params = state.focus_params or {}
            direction = params.get("direction") or (0.0, 0.0, -1.0)
            up = params.get("up") or (0.0, 1.0, 0.0)
            self._orient_camera(camera, direction, up)
            camera.frame_points(state.focus_points)
            # clip away whatever sits in front of the focused sphere, the way
            # Mol* sees into the cavity around a buried component. The sphere
            # is the component's, scaled/overridden by the focus parameters
            from ...framing import enclosing_sphere

            center, radius = enclosing_sphere(state.focus_points)
            if params.get("radius") is not None:
                radius = params["radius"] * _WORLD_SCALE
            else:
                radius = (
                    radius * params.get("radius_factor", 1.0)
                    + params.get("radius_extent", 0.0) * _WORLD_SCALE
                )
            camera.clip_to_sphere(center, radius)

    @staticmethod
    def _orient_camera(camera, forward, up, position=None) -> None:
        """Point the camera along ``forward`` with ``up`` roughly upward."""
        from mathutils import Matrix, Vector

        forward = Vector(forward).normalized()
        right = forward.cross(Vector(up)).normalized()
        if right.length == 0:
            # up parallel to the view direction; pick any perpendicular axis
            right = forward.orthogonal().normalized()
        true_up = right.cross(forward)

        rotation = Matrix((right, true_up, -forward)).transposed().to_4x4()
        if position is None:
            position = camera.camera.matrix_world.translation
        camera.camera.matrix_world = Matrix.Translation(position) @ rotation

    def _apply_camera(self, camera, params: dict) -> None:
        from mathutils import Vector

        position = Vector(params["position"]) * _WORLD_SCALE
        target = Vector(params["target"]) * _WORLD_SCALE
        up = Vector(params.get("up", (0.0, 1.0, 0.0)))
        self._orient_camera(camera, target - position, up, position=position)
        if params.get("near") is not None:
            camera.clip_start = max(params["near"] * _WORLD_SCALE, 1e-4)
