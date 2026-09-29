"""B-DNA (1BNA) to Z-DNA (2DCG), chain A.

After ProteinMotion's ``examples/dna_morph.py``. The C1' contact-guided morph
needs the planned ``Morph To Attribute`` node; the two structures are cut
between instead.
"""

from _common import hide, make_canvas, paint, render, show, spin, title, visibility
import molecularnodes as mn

canvas = make_canvas()
source = mn.Molecule.fetch("1BNA")
target = mn.Molecule.fetch("2DCG")
paint(source, "nucleic", "#75b8d9")
paint(target, "nucleic", "#df92a0")
source.add_style("cartoon", selection="chainID A", base_shape="Cylinder")
target.add_style("cartoon", selection="chainID A", base_shape="Cylinder")
target.add_style("ball_and_stick", selection="nucleic and chainID A")
sticks = target.styles["Style Ball and Stick"]
hide(sticks)
target.object.hide_render = True
canvas.look_at(source, viewpoint=(75, 0, -15), margin=0.05)
caption = title(source, "1BNA A → 2DCG A", position=(0.055, 0.89), size=27)
heading = title(source, "DNA morph · C1′ anchors", position=(0.055, 0.10), size=40)

with canvas.timeline() as t:
    show(t, caption, heading)
    t.wait(1)
    # PM: BackboneMorph(source, target, match=match, residue_delay=0.25, run_time=5)
    t.play(spin(source, 60), run_time=2.5)
    t.set((source.object, "hide_render"), True)
    t.set((target.object, "hide_render"), False)
    t.play(canvas.camera.look_at(target, margin=0.2), run_time=1.5)
    t.set(visibility(sticks), True)
    t.play(canvas.camera.orbit(20), run_time=2)
    t.wait(0.5)

render(canvas, t, "dna_morph")
