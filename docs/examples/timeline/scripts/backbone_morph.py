"""Calmodulin to troponin C.

After ProteinMotion's ``examples/backbone_morph.py``. The contact-guided
``BackboneMorph`` needs the planned ``Morph To Attribute`` node (a second set of
positions stored as an attribute with a keyable factor, staggered per residue).
Until then the two structures cross-fade.
"""

from _common import GOLD, TEAL, make_canvas, render, spin
import molecularnodes as mn

canvas = make_canvas()
source = mn.Molecule.fetch("1CLL")
target = mn.Molecule.fetch("1NCX")
source.set_color(TEAL)
target.set_color(GOLD)
source.add_style("cartoon", selection="chainID A")
target.add_style("cartoon", selection="chainID A")
canvas.look_at(source, viewpoint=(75, 0, 25), margin=0.25)

with canvas.timeline() as t:
    t.set_visible(target, False)
    t.wait(1)
    # PM: BackboneMorph(source, target, match=..., residue_delay=0.025, run_time=6)
    # t.play(source.morph(target, match, stagger=0.025), run_time=6)
    t.play(spin(source, 40), 3)
    t.play(
        [
            t.fade(source, 0.0),
            t.fade(target, 1.0),
            canvas.camera.look_at(target.get_view("chainID A"), margin=0.25),
        ],
        1.5,
    )
    t.play(spin(target, 100), 3)
    t.wait(0.5)

render(canvas, t, "backbone_morph")
