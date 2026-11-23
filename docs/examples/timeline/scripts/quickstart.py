"""A first film: a turn, a live region highlight, focus, and ball-and-stick.

After ProteinMotion's ``examples/quickstart.py`` (ubiquitin, 1UBQ).
"""

from _common import GOLD, make_canvas, render, switch
import molecularnodes as mn

canvas = make_canvas()
protein = mn.Molecule.fetch("1UBQ")
protein.set_color(GOLD, "resid 23-34")
protein.add_style("cartoon", selection="protein").add_style(
    "ball_and_stick", selection="protein"
)
cartoon, sticks = protein.styles[0], protein.styles[1]
helix = protein.get_view("resid 23-34")
canvas.look_at(protein, viewpoint=(70, 0, 20), margin=0.1)

with canvas.timeline() as t:
    t.set_visible(protein, False)
    # PM: FadeIn(protein) while the camera turns
    t.play(t.fade(cartoon, 1.0), canvas.camera.orbit(70), run_time=2)
    # PM: helix.highlight(style="box") - no highlight style yet, the helix is
    # coloured gold above instead.
    t.play(canvas.camera.look_at(helix, margin=0.3), run_time=1.5)
    t.play(canvas.camera.orbit(30), run_time=1.5)
    t.play(canvas.camera.look_at(protein, margin=0.1), run_time=1.2)
    switch(t, [cartoon, sticks], [sticks])
    t.wait(0.8)

render(canvas, t, "quickstart")
