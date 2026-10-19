"""Calmodulin 1CLL: select a helix, fade its context, then pull lens focus.

After ProteinMotion's ``examples/eevee_focus.py``. Depth of field is Blender's
own, with the focus empty the timeline's ``focus`` clip places.
"""

from _common import (
    BLUE,
    GOLD,
    callout,
    context_material,
    make_canvas,
    paint,
    render,
    show,
    title,
)
import molecularnodes as mn

canvas = make_canvas()
protein = mn.Molecule.fetch("1CLL")
paint(protein, None, BLUE)
paint(protein, "resid 5-19", GOLD)
helix = "protein and resid 5-19"
context = context_material(0.0)
protein.add_style("cartoon", selection=helix)
protein.add_style("cartoon", selection=f"protein and not ({helix})", material=context)
far_region = protein.get_view("resid 82-92 and name CA")
canvas.look_at(protein, viewpoint=(70, 0, 50), margin=0.12)
label = title(protein, "Calmodulin · 1CLL", position=(0.055, 0.1), size=30)
note = callout(
    protein, helix, "Alpha helix|Chain A · residues 5–19", offset=(220, -90), size=34
)

with canvas.timeline() as t:
    t.play(
        canvas.camera.focus(protein.get_view(f"{helix} and name CA"), fstop=1.8),
        run_time=0,
    )
    show(t, label, note)  # PM: Write(label) - no write-on yet
    # PM: SetOpacity(context, 0.06) - the context is its own style with a
    # Transparent material whose transparency socket is tweened
    t.play(t.tween(context.node.i.transparency, 0.94), run_time=1)
    t.play(canvas.camera.orbit(16), run_time=2)
    t.wait(0.4)
    t.play(
        canvas.camera.focus(far_region, fstop=1.2),
        t.tween(context.node.i.transparency, 0.65),
        run_time=1.4,
    )
    t.wait(0.4)
    t.play(
        canvas.camera.focus(protein.get_view(f"{helix} and name CA"), fstop=1.8),
        t.tween(context.node.i.transparency, 0.94),
        run_time=1.2,
    )
    t.wait(0.8)

render(canvas, t, "eevee_focus")
