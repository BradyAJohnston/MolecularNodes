"""Ubiquitin NMR conformers (2K39): playback with camera focus on regions.

After ProteinMotion's ``examples/nmr_regions.py``. The 3D region highlights
(sphere, box, atom halos) have no style node yet; the regions are painted.
"""

from _common import (
    CYAN,
    GOLD,
    align_frames,
    hide,
    make_canvas,
    paint,
    render,
    show_only,
)
import molecularnodes as mn

canvas = make_canvas()
protein = mn.Molecule.fetch("2K39")
align_frames(protein, "name CA and resid 1-70")
paint(protein, "resid 23-34", GOLD)
paint(protein, "resid 71-76", CYAN)
protein.add_style("cartoon").add_style("ball_and_stick")
cartoon, sticks = list(protein.styles)
hide(sticks)
helix = protein.get_view("resid 23-34")
tail = protein.get_view("resid 71-76")
canvas.look_at(protein, viewpoint=(70, 0, 40), margin=0.15)

with canvas.timeline() as t:
    t.wait(0.5)
    # PM: PlayTrajectory(p, start=0, end=20, state_easing=smooth). The entity
    # interpolates between models when its `interpolate` property is on.
    t.play(protein.play(0, 20), run_time=4, easing="smooth")
    # PM: FadeIn(sphere), FadeIn(box) around the helix
    t.play(canvas.camera.look_at(helix, margin=0.3), run_time=1.8)
    t.play(protein.play(20, 35), run_time=4.5, easing="smooth")
    t.play(canvas.camera.look_at(protein, margin=0.15), run_time=1.8)
    # PM: FadeIn(halo), FadeIn(tail_box) around the tail
    t.play(canvas.camera.look_at(tail, margin=0.3), run_time=1.8)
    show_only(t, [cartoon, sticks], [sticks])
    t.play(protein.play(35, 45), run_time=5, easing="smooth")
    t.play(canvas.camera.look_at(protein, margin=0.15), run_time=1.8)
    t.wait(0.5)

render(canvas, t, "nmr_regions")
