"""Map 1UBQ B factors to residue colour and cartoon thickness.

After ProteinMotion's ``examples/numerical_properties.py``.
"""

from _common import INK, make_canvas, render, show, spin, switch, title
import molecularnodes as mn
from molecularnodes.nodes import geometry as mg

canvas = make_canvas()
protein = mn.Molecule.fetch("1UBQ", hide="water")
mat = mn.material.Default().material

with protein.tree.reset() as (atoms, join):
    atoms = (
        atoms
        >> mg.SetColor(
            color=mg.ColorMatplotlib(mg.FieldRemap(mg.BFactor()), category="Diverging")
        )
    ).o.atoms

    cartoon = mg.StyleCartoon(atoms, material=mat)
    sticks = mg.StyleBallAndStick(atoms, material=mat)
    surface = mg.StyleSurface(atoms, material=mat)

    _ = list(x >> join for x in (cartoon, sticks, surface))


canvas.look_at(protein, viewpoint=(70, 0, 30), margin=0.2)
heading = title(protein, "Residue properties", size=42)
detail = title(
    protein, "Ubiquitin · Cα B factors", position=(0.06, 0.13), size=25, color=INK
)
# PM: ColorLegend(scale, title="B factor", unit="Å²") - no legend annotation yet.

with canvas.timeline() as t:
    t.set_visible([sticks, surface], False)
    show(t, heading, detail)
    # PM: ColorByProperty(..., thickness=(0.6, 1.8), residue_delay=0.012). Colour
    # by field is the node above; per-residue thickness by field and the
    # staggered reveal need the planned "cartoon thickness by field" and
    # "Animate Stagger" nodes, so the whole loop radius thickens instead.
    t.play((t.tween(cartoon.i.loop_radius, 0.5), spin(protein, 15)), 2.5)
    t.wait(0.5)
    switch(t, [cartoon, sticks, surface], [sticks])
    t.wait(0.5)
    t.play(spin(protein, 17), 1)
    switch(t, [cartoon, sticks, surface], [surface])
    t.wait(0.5)
    t.play(spin(protein, 15), 1.5)

render(canvas, t, "numerical_properties")
