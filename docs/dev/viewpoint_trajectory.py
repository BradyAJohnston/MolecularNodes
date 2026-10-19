"""
Render Orient to Best View over the AdK open-to-closed transition (MDAnalysisData).

    uv run docs/dev/viewpoint_trajectory.py OUT_DIR [VARIANT ...]
    uv run docs/dev/viewpoint_trajectory.py OUT_DIR combine

Renders the same playback with a fixed front view and with Orient to Best View
under the settings in VARIANTS: jumping to each frame's best view, eased with a
half life, held near the current view with stickiness, or both, and a sweep of
stickiness levels. The target is the active-site residues, which the LID and NMP
domains close over. The protein is drawn as a cartoon, as a surface, or as a
cartoon with the side chains near the site, and that drawing is what can block
the view.

Naming variants renders only those and keeps earlier results in OUT_DIR for the
rest; ``combine`` only re-tiles the MP4s already there. Writes one MP4 per
variant, a labelled grid per entry in GRIDS, and results.json with, per
frame, how many target atoms have a clear ray to the camera and how far the view
turned since the previous frame.
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())
import bpy  # noqa: E402
import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402
from MDAnalysisData import datasets  # noqa: E402
from nodebpy.nodes import geometry as g  # noqa: E402
import molecularnodes as mn  # noqa: E402
from molecularnodes import viewpoint as vp  # noqa: E402
from molecularnodes.nodes.geometry import (  # noqa: E402
    IsSideChain,
    OrientToBestView,
    SelectProximity,
    SetColor,
    StyleBallAndStick,
    StyleCartoon,
    StyleSticks,
    StyleSurface,
)
from molecularnodes.scene.base import _WorldBVH  # noqa: E402

# P-loop lysine, the arginines that bind the phosphates, and Lys200
ACTIVE_SITE = "resid 13 36 88 123 156 167 200"
VARIANTS = {
    "fixed": dict(orient=False),
    "jump": dict(orient=True, smooth=False),
    "smooth": dict(orient=True, smooth=True, half_life=0.5),
    "smooth_1s": dict(orient=True, smooth=True, half_life=1.0),
    "sticky": dict(orient=True, stickiness=0.25),
    "sticky_smooth": dict(orient=True, smooth=True, half_life=0.5, stickiness=0.25),
    **{f"sticky_{s}": dict(orient=True, stickiness=s) for s in (0.05, 0.1, 0.5, 1.0)},
    **{
        f"{style}_{name}": dict(style=style, **options)
        for style in ("surface", "sidechains")
        for name, options in (
            ("fixed", dict(orient=False)),
            (
                "sticky_smooth",
                dict(orient=True, smooth=True, half_life=0.5, stickiness=0.25),
            ),
        )
    },
}
LABELS = {
    "fixed": "fixed front view",
    "jump": "orient",
    "smooth": "orient, half life 0.5 s",
    "smooth_1s": "orient, half life 1.0 s",
    "sticky": "orient, stickiness 0.25",
    "sticky_smooth": "orient, stickiness 0.25, half life 0.5 s",
    **{f"sticky_{s}": f"orient, stickiness {s}" for s in (0.05, 0.1, 0.5, 1.0)},
    **{
        f"{style}_{name}": f"{label}, {view}"
        for style, label in (("surface", "surface"), ("sidechains", "side chains"))
        for name, view in (
            ("fixed", "fixed"),
            ("sticky_smooth", "stickiness 0.25, half life 0.5 s"),
        )
    },
}
# each grid tiles its variants, in order, into OUT_DIR/<grid>.mp4
GRIDS = {
    "side_by_side": ["fixed", "jump", "smooth", "smooth_1s", "sticky", "sticky_smooth"],
    "stickiness": [
        "jump",
        "sticky_0.05",
        "sticky_0.1",
        "sticky",
        "sticky_0.5",
        "sticky_1.0",
    ],
    "styles": [
        "fixed",
        "surface_fixed",
        "sidechains_fixed",
        "sticky_smooth",
        "surface_sticky_smooth",
        "sidechains_sticky_smooth",
    ],
}
NEAR = 5.0  # Å from the site for side chains to show
FONT = "/usr/share/fonts/TTF/DejaVuSansMono.ttf"
FPS = 24
SUBFRAMES = 1


def turn(a, b):
    """Degrees of the best-fit rotation taking points a onto points b (Kabsch)."""
    a, b = a - a.mean(axis=0), b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1.0, 1.0, d]) @ vt
    return float(np.degrees(np.arccos(np.clip((np.trace(r) - 1) / 2, -1, 1))))


def build(mol, orient, smooth=False, half_life=0.5, stickiness=0.0, style="cartoon"):
    material = mn.material.Default().material
    with mol.tree.reset() as (atoms, join):
        is_site = g.NamedAttribute.boolean("is_site")
        grey = atoms >> SetColor(color=(0.6, 0.6, 0.65, 1.0))
        if style == "surface":
            # the site is drawn on its own, so leave it out of the surface
            protein = grey >> StyleSurface(selection=~is_site, material=material)
        else:
            protein = grey >> StyleCartoon(material=material)
        site = (
            atoms
            >> SetColor(color=(0.95, 0.55, 0.15, 1.0))
            >> StyleBallAndStick(selection=is_site, material=material)
        )
        parts = [protein, site]
        if style == "sidechains":
            near = SelectProximity(
                target_atoms=g.SeparateGeometry(atoms, selection=is_site).o.selection,
                expand=True,
                distance_a=NEAR,
            ).o.selection
            side = near & ~is_site & IsSideChain().o.selection
            parts.append(atoms >> StyleSticks(selection=side, material=material))
        styled = g.JoinGeometry(parts)
        moved = atoms
        if orient:
            node = OrientToBestView(
                geometry=styled,
                atoms=atoms,
                selection=g.NamedAttribute.boolean("is_site"),
                smooth=smooth,
                half_life=half_life,
                stickiness=stickiness,
            )
            styled = node.o.geometry
            moved = atoms >> g.TransformGeometry(
                transform=node.o.transform, mode="Matrix"
            )
        # loose atom vertices carry where the atoms ended up and render as nothing
        probe = moved >> g.StoreNamedAttribute.point.vector(
            name="moved", value=g.Position().o.position
        )
        g.JoinGeometry([probe, styled]) >> join


def combine(out: Path, grid: str, columns: int = 3):
    """Tile the rendered variants of a grid into one labelled MP4."""
    rendered = [v for v in GRIDS[grid] if (out / f"{v}.mp4").exists()]
    labelled = ";".join(
        f"[{i}:v]drawtext=fontfile={FONT}:text='{LABELS[v]}':x=10:y=10"
        f":fontsize=16:fontcolor=white:box=1:boxcolor=black@0.5:boxborderw=4[v{i}]"
        for i, v in enumerate(rendered)
    )
    layout = "|".join(
        "_".join(
            "+".join(f"{axis}0" for _ in range(k)) or "0"
            for axis, k in (("w", i % columns), ("h", i // columns))
        )
        for i in range(len(rendered))
    )
    stack = "".join(f"[v{i}]" for i in range(len(rendered)))
    subprocess.run(
        ["ffmpeg", "-y", "-loglevel", "error"]
        + [a for v in rendered for a in ("-i", str(out / f"{v}.mp4"))]
        + [
            "-filter_complex",
            f"{labelled};{stack}xstack=inputs={len(rendered)}"
            f":layout={layout}:fill=black",
        ]
        + ["-pix_fmt", "yuv420p", "-c:v", "libx264", "-crf", "26"]
        + [str(out / f"{grid}.mp4")],
        check=True,
    )


def main(out: Path, names: list[str]):
    if names == ["combine"]:
        for grid in GRIDS:
            combine(out, grid)
        return
    out.mkdir(parents=True, exist_ok=True)
    data = datasets.fetch_adk_transitions_DIMS()
    canvas = mn.Canvas(engine="EEVEE", resolution=(480, 360))
    canvas.fps = FPS
    saved = out / "results.json"
    results = json.loads(saved.read_text()) if saved.exists() else {}
    for name in names or VARIANTS:
        options = VARIANTS[name]
        canvas.clear()
        u = mda.Universe(data.topology, data.trajectories[0])
        mol = mn.Molecule(u)
        mol.dssp.init()
        mol.subframes, mol.interpolate = SUBFRAMES, True
        site = np.zeros(len(mol), dtype=bool)
        site[u.select_atoms(ACTIVE_SITE).indices] = True
        mol.store_named_attribute(site, "is_site")
        build(mol, **options)
        n_frames = u.trajectory.n_frames * (SUBFRAMES + 1) - 1
        canvas.frame_range = (0, n_frames - 1)
        canvas.frame = 0
        n = len(mol)
        # the oriented variants turn about the site, so frame a sphere about it
        # that holds the protein in any orientation, and the same for every variant
        first = mol.named_attribute("moved", evaluate=True)[:n]
        pivot = first[site].mean(axis=0)
        radius = np.linalg.norm(first - pivot, axis=1).max()
        canvas.look_at(pivot + vp.fibonacci_sphere(200) * radius, margin=0.02)
        clear, turns, previous = [], [], first
        with canvas.record(str(out / f"{name}.mp4"), fps=FPS) as movie:
            # stepping one frame at a time is what advances the smoothing
            for frame in range(n_frames):
                canvas.frame = frame
                moved = mol.named_attribute("moved", evaluate=True)[:n]
                forward = -np.array(canvas.camera.camera.matrix_world.col[2][:3])
                tree = _WorldBVH(mol.object, bpy.context.evaluated_depsgraph_get())
                clear.append(int(vp.ray_scorer([tree], moved[site])(forward[None])[0]))
                turns.append(turn(previous, moved))
                previous = moved
                movie.render()
        results[name] = dict(clear=clear, turn=turns)
        print(name, "mean clear", np.mean(clear), "of", int(site.sum()), flush=True)
    results["atoms"] = int(site.sum())
    results["fps"] = FPS
    (out / "results.json").write_text(json.dumps(results))
    for grid in GRIDS:
        combine(out, grid)


if __name__ == "__main__":
    main(Path(sys.argv[1]), sys.argv[2:])
