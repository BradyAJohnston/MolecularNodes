"""B-DNA (1BNA) to Z-DNA (2DCG), chain A.

After ProteinMotion's ``examples/dna_morph.py``. The C1' contact-guided morph
needs the planned ``Morph To Attribute`` node; the two structures are
cross-faded instead.
"""

from _common import make_canvas, render, show, spin, title
import molecularnodes as mn

canvas = make_canvas()
source = mn.Molecule.fetch("1BNA")
target = mn.Molecule.fetch("2DCG")
source.set_color("#75b8d9", "nucleic")
target.set_color("#df92a0", "nucleic")
source.add_style("cartoon", selection="chainID A", base_shape="Cylinder")
target.add_style("cartoon", selection="chainID A", base_shape="Cylinder")
target.add_style("ball_and_stick", selection="nucleic and chainID A")
sticks = target.styles["Style Ball and Stick"]
canvas.look_at(source, viewpoint=(75, 0, -15), margin=0.05)
caption = title(source, "1BNA A → 2DCG A", position=(0.055, 0.89), size=27)
heading = title(source, "DNA morph · C1′ anchors", position=(0.055, 0.10), size=40)

with canvas.timeline() as t:
    t.set_visible(target, False)
    show(t, caption, heading)
    t.wait(1)
    # PM: BackboneMorph(source, target, match=match, residue_delay=0.25, run_time=5)
    t.play(spin(source, 60), run_time=2.5)
    t.play(
        t.fade(source, 0.0),
        t.fade(target.styles[0], 1.0),
        canvas.camera.look_at(target.get_view("chainID A"), margin=0.2),
        run_time=1.5,
    )
    t.play(t.fade(sticks, 1.0), run_time=0.5)
    t.play(canvas.camera.orbit(20), run_time=2)
    t.wait(0.5)

render(canvas, t, "dna_morph")
