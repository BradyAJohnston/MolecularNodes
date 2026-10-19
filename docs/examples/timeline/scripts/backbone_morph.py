"""Calmodulin to troponin C.

After ProteinMotion's ``examples/backbone_morph.py``. The contact-guided
``BackboneMorph`` needs the planned ``Morph To Attribute`` node (a second set of
positions stored as an attribute with a keyable factor, staggered per residue).
Until then the two structures are shown one after the other.
"""

from _common import GOLD, TEAL, make_canvas, paint, render, spin
import molecularnodes as mn

canvas = make_canvas()
source = mn.Molecule.fetch("1CLL")
target = mn.Molecule.fetch("1NCX")
paint(source, None, TEAL)
paint(target, None, GOLD)
source.add_style("cartoon", selection="chainID A")
target.add_style("cartoon", selection="chainID A")
# every entity is one object; a cut is the target style switched on and the
# source switched off. Both styles are on selections, so key the objects'
# render visibility instead of the style nodes.
canvas.look_at(source, viewpoint=(75, 0, 25), margin=0.25)
target.object.hide_render = True

with canvas.timeline() as t:
    t.wait(1)
    # PM: BackboneMorph(source, target, match=..., residue_delay=0.025, run_time=6)
    # t.play(source.morph(target, match, stagger=0.025), run_time=6)
    t.play(spin(source, 40), run_time=3)
    t.set((source.object, "hide_render"), True)
    t.set((target.object, "hide_render"), False)
    t.play(canvas.camera.look_at(target, margin=0.25), run_time=1.5)
    t.play(spin(target, 100), run_time=3)
    t.wait(0.5)

render(canvas, t, "backbone_morph")
