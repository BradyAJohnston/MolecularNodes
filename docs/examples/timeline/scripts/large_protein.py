"""Actual GroEL/GroES coordinates: 21 chains, coloured by chain.

After ProteinMotion's ``examples/large_protein.py``.
"""

from _common import hide, make_canvas, render, show_only, spin
import molecularnodes as mn
from molecularnodes.nodes import geometry as mg

canvas = make_canvas()
complex_ = mn.Molecule.fetch("1AON")
complex_.add_style("cartoon", color=lambda: mg.ColorAttributeRandom(name="chain_id"))
complex_.add_style(
    "ball_and_stick", color=lambda: mg.ColorAttributeRandom(name="chain_id")
)
cartoon, sticks = complex_.styles[0], complex_.styles[1]
hide(sticks)
# PM orients the chaperonin axis vertically from an SVD of the Cα positions;
# the side viewpoint does the same job for this assembly.
canvas.look_at(complex_, viewpoint=(80, 0, 17), margin=0.05)
# PM: camera.depth_cue = 0.45 - no depth cue (mist) shortcut yet.

with canvas.timeline() as t:
    t.wait(0.5)
    t.play(spin(complex_, 180), run_time=4, easing="linear")
    show_only(t, [cartoon, sticks], [sticks])
    t.play(spin(complex_, 90), run_time=3, easing="linear")
    t.wait(0.5)

render(canvas, t, "large_protein_complex")

# a single subunit, on its own
canvas.clear()
subunit = mn.Molecule.fetch("1AON")
subunit.add_style("cartoon", selection="chainID A")
canvas.look_at(subunit, viewpoint=(80, 0, 17), margin=0.05)
with canvas.timeline() as t:
    t.play(spin(subunit, 360), run_time=6, easing="linear")

render(canvas, t, "large_protein_subunit")
