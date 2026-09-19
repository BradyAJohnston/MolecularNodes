"""
Demo storyboard for the Canvas timeline prototype.

Loads a trajectory, then plays: a framing move, an orbit while a ribbon
thickens, a focus pull with depth of field onto a selection, a trajectory
playback with a dolly-in, and a final hold. Renders to MP4 and a contact sheet.
"""

import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())

import bpy  # noqa: E402
from nodebpy.nodes.geometry import NamedAttribute  # noqa: E402
import molecularnodes as mn  # noqa: E402
from molecularnodes.nodes import geometry as mg  # noqa: E402

out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
out_dir.mkdir(parents=True, exist_ok=True)
data_dir = Path("tests/data")

canvas = mn.Canvas(engine="EEVEE", resolution=(960, 540))
canvas.engine.samples = 32
canvas.background = (0.08, 0.08, 0.1, 1.0)

traj = mn.Molecule.load(
    data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
)
with traj.tree as tree:
    ribbon = mg.StyleRibbon(material=mn.material.AmbientOcclusion().material)
    tree.atoms >> mg.SetColor(color=mg.ColorRainbow()) >> ribbon >> tree.join
    sticks = mg.StyleBallAndStick(
        selection=NamedAttribute.boolean("is_side_chain"),
        material=mn.material.Default().material,
    )
    tree.atoms >> mg.SetColor(color=mg.ColorCommon()) >> sticks >> tree.join

label = traj.annotations.add_label_2d(text="Storyboard demo", location=(0.05, 0.92))
label.text_size = 28

canvas.look_at(traj, viewpoint="front", margin=0.1)
ribbon.i.peptide_radius.socket.default_value = 0.6
site = traj.get_view("resid 40-60")

with canvas.timeline(fps=24) as t:
    # 0-1.5 s: settle into a wider framing from the side
    t.play(
        canvas.camera.look_at(traj, viewpoint=(70, 0, 35), margin=0.15), run_time=1.5
    )
    # 1.5-4.5 s: sweep round while the ribbon fattens
    t.play(
        canvas.camera.orbit(120, about=traj),
        t.tween(ribbon.i.peptide_radius, 2.0),
        run_time=3,
    )
    t.wait(0.25)
    # pull focus onto the site with a shallow depth of field and push in
    t.play(
        canvas.camera.focus(site, fstop=0.8),
        canvas.camera.dolly(1.5),
        run_time=1.5,
    )
    # play the trajectory while the label shrinks away
    t.play(
        t.frames(traj, 0, 4),
        t.tween((label, "text_size"), 8),
        run_time=2,
    )
    t.wait(0.5)

print(t)
t.render(out_dir / "storyboard.mp4")

# a contact sheet of a few frames for the write-up
from PIL import Image  # noqa: E402

frames = [t.start + int(f) for f in (0, 36, 72, 108, 144, 180, 216)]
frames = [f for f in frames if f <= t.frame_end]
tiles = []
for frame in frames:
    canvas.snapshot(out_dir / f"frame_{frame:04d}.png", frame=frame, render_scale=50)
    tiles.append(Image.open(out_dir / f"frame_{frame:04d}.png"))
w, h = tiles[0].size
sheet = Image.new("RGB", (w * 4, h * ((len(tiles) + 3) // 4)), (20, 20, 24))
for i, tile in enumerate(tiles):
    sheet.paste(tile, ((i % 4) * w, (i // 4) * h))
sheet.save(out_dir / "contact_sheet.png")
bpy.ops.wm.save_as_mainfile(filepath=str(out_dir / "storyboard.blend"))
print("done", out_dir)
