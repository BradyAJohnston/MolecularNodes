"""Yeast tRNA: base shapes, a loop callout, B-factor colours and a surface.

After ProteinMotion's ``examples/rna_styles.py`` (1EHZ).
"""

from _common import (
    GOLD,
    callout,
    context_material,
    hide,
    make_canvas,
    paint,
    render,
    show,
    show_only,
    title,
    visibility,
)
import molecularnodes as mn
from molecularnodes.nodes import geometry as mg

canvas = make_canvas()
rna = mn.Molecule.fetch("1EHZ")
loop = "nucleic and resid 32-38"
paint(rna, loop, GOLD)
context = context_material(0.0)
# PM's BaseStyle "rings" / "slabs" are the cartoon's base shapes; each is a
# style node here and the timeline cuts between them
rna.add_style("cartoon", selection=loop, base_shape="Rectangle")
rna.add_style(
    "cartoon",
    selection=f"nucleic and not ({loop})",
    material=context,
    base_shape="Rectangle",
)
rna.add_style("cartoon", selection="nucleic", base_shape="Cylinder")
rna.add_style(
    "cartoon",
    selection="nucleic",
    color=lambda: mg.ColorAttributeMap(name="b_factor", min=20, max=90),
)
rna.add_style(
    "surface",
    selection="nucleic",
    color=lambda: mg.ColorAttributeMap(name="b_factor", min=20, max=90),
)
slabs, faint, rings, by_bfactor, surface = (
    rna.styles["Style Cartoon"],
    rna.styles["Style Cartoon.001"],
    rna.styles["Style Cartoon.002"],
    rna.styles["Style Cartoon.003"],
    rna.styles["Style Surface"],
)
hide(rings, by_bfactor, surface)
canvas.look_at(rna, viewpoint=(70, 0, -30), margin=0.2)
heading = title(rna, "Transfer RNA · 1EHZ", position=(0.055, 0.10), size=30)
note = callout(rna, loop, "Anticodon loop|Residues 32–38", offset=(200, -60), size=34)
# PM: ColorLegend(scale, title="Deposited B factor") - no legend annotation yet.

with canvas.timeline() as t:
    t.play(canvas.camera.focus(rna.get_view(loop), fstop=5.6), run_time=0)
    show(t, heading)
    t.play(canvas.camera.orbit(20), run_time=3)
    # PM: Colorize(loop, gold, residue_delay=0.08), SetOpacity(context, 0.18)
    t.play(t.tween(context.node.i.transparency, 0.82), run_time=1.5)
    show(t, note)
    t.play(canvas.camera.orbit(9), run_time=1.5)
    # PM: BaseStyle(rna, "slabs") -> here a cut to base rings, showing every
    # residue with the same style again
    t.set(visibility(rings), True)
    t.play(canvas.camera.orbit(9), run_time=1.5)
    show(t, note, visible=False)
    t.play(t.tween(context.node.i.transparency, 0.0), run_time=1)
    show_only(t, [slabs, faint, rings, by_bfactor, surface], [by_bfactor])
    t.wait(1.5)
    show_only(t, [slabs, faint, rings, by_bfactor, surface], [surface])
    t.play(canvas.camera.orbit(9), run_time=1.5)
    t.play(canvas.camera.orbit(17), run_time=2.5)
    t.wait(0.5)

render(canvas, t, "rna_styles")
