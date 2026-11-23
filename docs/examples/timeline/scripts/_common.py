"""
Shared setup for the timeline examples.

Every example is a stand-alone script: ``uv run python docs/examples/timeline/scripts/<name>.py
[out_dir]``. Set ``MN_EXAMPLES_QUICK=1`` to render small and fast (for checking
that a script runs), and ``MN_EXAMPLES_ENGINE=CYCLES`` on a machine without a
GPU.
"""

import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())

import bpy  # noqa: E402
import numpy as np  # noqa: E402
import molecularnodes as mn  # noqa: E402

QUICK = os.environ.get("MN_EXAMPLES_QUICK", "") not in ("", "0")
ENGINE = os.environ.get("MN_EXAMPLES_ENGINE", "EEVEE").upper()

# the palette ProteinMotion's examples use, so the scenes read the same
BLUE, GOLD, TEAL = "#759dbb", "#efb766", "#64d6cb"
CYAN, CORAL, INK, WHITE = "#50e0d0", "#f06478", "#94a6ba", "#edf1f6"
BACKGROUND = (0.03, 0.06, 0.10, 1.0)


def output_dir() -> Path:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("examples_out")
    out.mkdir(parents=True, exist_ok=True)
    return out


def make_canvas(fps: int = 24, resolution=(960, 540)) -> mn.Canvas:
    """A canvas with the dark studio look of the ProteinMotion examples."""
    if QUICK:
        resolution, fps = (320, 180), 6
    if ENGINE == "CYCLES":
        engine = mn.scene.Cycles(samples=4 if QUICK else 64, device="CPU")
    else:
        engine = mn.scene.EEVEE(samples=4 if QUICK else 32)
    canvas = mn.Canvas(engine, resolution=resolution)
    canvas.compositor.device = "CPU"
    canvas.fps = fps
    canvas.background = BACKGROUND
    clear_annotation_overlay()
    return canvas


def clear_annotation_overlay() -> None:
    """
    Blank the compositor's annotation overlay image.

    The render handler only redraws that image when some entity has
    annotations, so after ``canvas.clear()`` the last scene's labels would be
    composited over every later render (a Molecular Nodes bug, seen
    2026-09-21). Scripts that render several scenes in one process call this
    between them.
    """
    from molecularnodes.handlers import annotations_image

    image = bpy.data.images.get(annotations_image)
    if image is not None:
        image.pixels.foreach_set(np.zeros(len(image.pixels), dtype=np.float32))
        image.update()


def render(canvas: mn.Canvas, timeline, name: str, frames: int = 4) -> Path:
    """Render the storyboard to ``<out>/<name>.mp4`` plus a contact sheet."""
    out = output_dir()
    path = out / f"{name}.mp4"
    timeline.render(path)
    picks = np.linspace(timeline.start, timeline.frame_end, frames).astype(int)
    images = []
    for frame in picks:
        image = canvas.snapshot(out / f"{name}_{frame:04d}.png", frame=int(frame))
        images.append(image)
    try:
        from PIL import Image

        tiles = [Image.open(out / f"{name}_{f:04d}.png") for f in picks]
        w, h = tiles[0].size
        sheet = Image.new("RGB", (w * len(tiles), h))
        for i, tile in enumerate(tiles):
            sheet.paste(tile.convert("RGB"), (i * w, 0))
        sheet.save(out / f"{name}_sheet.png")
    except ImportError:
        pass
    print(f"wrote {path}")
    print(timeline)
    return path


def rgba(hex_color: str, alpha: float = 1.0) -> tuple[float, float, float, float]:
    """A CSS hex colour as linear RGBA."""
    h = hex_color.lstrip("#")
    srgb = np.array([int(h[i : i + 2], 16) / 255 for i in (0, 2, 4)])
    linear = np.where(srgb <= 0.04045, srgb / 12.92, ((srgb + 0.055) / 1.055) ** 2.4)
    return (*map(float, linear), alpha)


def switch(timeline, styles, shown, run_time: float = 0.5) -> None:
    """
    Cross-fade from one representation to another: ProteinMotion's
    ``Representation(p, name)``. Every style in ``styles`` fades to 1 if it is
    in ``shown`` and to 0 otherwise.
    """
    # nodebpy overloads ``==`` on nodes to build Compare nodes, so test identity
    clips = [
        timeline.fade(style, 1.0 if any(style is s for s in shown) else 0.0)
        for style in styles
    ]
    timeline.play(*clips, run_time=run_time)


def align_frames(mol: mn.Molecule, selection: str = "name CA") -> None:
    """Superpose every frame on the first over ``selection`` (in memory)."""
    from MDAnalysis.analysis import align

    u = mol.universe
    align.AlignTraj(u, u, select=selection, in_memory=True).run()


def turntable(mol: mn.Molecule) -> bpy.types.Object:
    """
    An empty at the centre of the molecule's rendered geometry with the
    object parented to it, so a rotation of the empty turns the molecule
    about its own centre rather than its data origin. Made on first use.
    """
    from molecularnodes.blender import utils

    obj = mol.object
    if obj.parent is not None and obj.parent.name.endswith(" turntable"):
        return obj.parent
    from mathutils import Matrix

    centre = Matrix.Translation(utils.evaluated_points(obj).mean(axis=0))
    empty = bpy.data.objects.new(f"{obj.name} turntable", None)
    empty.empty_display_type = "PLAIN_AXES"
    empty.matrix_world = centre
    bpy.context.scene.collection.objects.link(empty)
    obj.parent = empty
    # the parent inverse cancels the empty's offset, so the molecule stays
    # put; a fresh empty's matrix_world is stale until the depsgraph runs,
    # hence the matrix is built rather than read back
    obj.matrix_parent_inverse = centre.inverted()
    return empty


def spin(mol: mn.Molecule, degrees: float, axis: str = "z", run_time=None, easing=None):
    """
    A clip turning the molecule about a world axis through its centre,
    ProteinMotion's ``Rotate(protein, angle)``. Successive spins accumulate.
    """
    import math
    from molecularnodes.scene.timeline import Channel, Tween

    empty = turntable(mol)
    index = "xyz".index(axis)
    channel = Channel(empty, "rotation_euler", index)
    end = channel.get() + math.radians(degrees)
    return Tween(pairs=[(channel, end)], run_time=run_time, easing=easing)


def design_pixels(size: float) -> int:
    """
    Text sizes in the examples are ProteinMotion's, in pixels of a 1080p
    frame; annotation text is in render pixels, so scale by the frame height.
    """
    # a label's glyphs render about 1.44 times its nominal size (measured on
    # a 540-line frame), so divide to match ProteinMotion's visual size
    return int(round(size * bpy.context.scene.render.resolution_y / 1080 / 1.44))


def title(mol: mn.Molecule, text: str, position=(0.06, 0.07), size=36, color=WHITE):
    """
    A 2D text annotation, placed as ProteinMotion places ``Text``: fractions
    of the frame from the top-left. Starts hidden; step ``visible`` to show it.
    """
    # ProteinMotion positions the top-left corner of the text; a 2D label is
    # placed by its bottom-left, so drop it by the text height (design pixels
    # of a 1080-line frame, multi-line text counted per line)
    lines = text.count("|") + 1
    y = 1.0 - position[1] - size * lines / 1080
    label = mol.annotations.add_label_2d(text=text, location=(position[0], y))
    label.text_size = design_pixels(size)
    label.text_color = rgba(color)
    label.text_align = "left"
    label.visible = False
    return label


def callout(
    mol: mn.Molecule, selection: str, text: str, offset=(180, 60), size=26, color=GOLD
):
    """
    A 3D label anchored at a selection's centre, with a short pointer and the
    text offset in pixels: ProteinMotion's ``region.callout``. Starts hidden.
    Annotation positions are in angstrom (the universe's frame), not world units.
    """
    centre = mol.universe.select_atoms(selection).positions.mean(axis=0)
    label = mol.annotations.add_label_3d(text=text, location=tuple(map(float, centre)))
    label.text_size = design_pixels(size)
    label.text_color = rgba(color)
    label.line_color = rgba(color)
    label.line_pointer_length = 6.0
    label.text_depth = False
    label.text_offset_x, label.text_offset_y = (design_pixels(v) for v in offset)
    label.visible = False
    return label


def show(timeline, *labels, visible: bool = True) -> None:
    """Step annotations visible (or hidden) at the cursor: a cut, not a write-on."""
    for label in labels:
        timeline.set((label, "visible"), visible)


def world_points(*mols: mn.Molecule) -> np.ndarray:
    """The rendered geometry of several entities, in world space, for framing."""
    from molecularnodes.blender import utils

    return np.concatenate([utils.evaluated_points(mol.object) for mol in mols])
