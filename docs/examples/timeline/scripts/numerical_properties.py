"""Map 1UBQ B factors to residue colour and cartoon thickness.

After ProteinMotion's ``examples/numerical_properties.py``.
"""

from _common import INK, hide, make_canvas, render, show, show_only, spin, title
import molecularnodes as mn
from molecularnodes.nodes import geometry as mg

canvas = make_canvas()
protein = mn.Molecule.fetch("1UBQ")
# the colour scale, min 0 to max 40 Å², is the ``Color Attribute Map`` node
by_bfactor = lambda: mg.ColorAttributeMap(name="b_factor", min=0, max=40)  # noqa: E731
protein.add_style("cartoon", color=by_bfactor)
protein.add_style("ball_and_stick", color=by_bfactor)
protein.add_style("surface", color=by_bfactor)
cartoon, sticks, surface = protein.styles[0], protein.styles[1], protein.styles[2]
hide(sticks, surface)
canvas.look_at(protein, viewpoint=(70, 0, 30), margin=0.2)
heading = title(protein, "Residue properties", size=42)
detail = title(
    protein, "Ubiquitin · Cα B factors", position=(0.06, 0.13), size=25, color=INK
)
# PM: ColorLegend(scale, title="B factor", unit="Å²") - no legend annotation yet.

with canvas.timeline() as t:
    show(t, heading, detail)
    # PM: ColorByProperty(..., thickness=(0.6, 1.8), residue_delay=0.012). Colour
    # by field is the node above; per-residue thickness by field and the
    # staggered reveal need the planned "cartoon thickness by field" and
    # "Animate Stagger" nodes, so the whole loop radius thickens instead.
    t.play(t.tween(cartoon.i.loop_radius, 0.5), spin(protein, 15), run_time=2.5)
    t.wait(0.5)
    show_only(t, [cartoon, sticks, surface], [sticks])
    t.play(spin(protein, 17), run_time=1)
    show_only(t, [cartoon, sticks, surface], [surface])
    t.play(spin(protein, 15), run_time=1.5)

render(canvas, t, "numerical_properties")
