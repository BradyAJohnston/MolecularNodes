"""Titles, live region callouts and residue labels on the 2K39 NMR ensemble.

After ProteinMotion's ``examples/labels_and_callouts.py``. Write / Unwrite
(vector write-on) needs a ``progress`` property on text annotations; labels
are cut in and out here.
"""

from _common import (
    CYAN,
    GOLD,
    align_frames,
    callout,
    make_canvas,
    render,
    show,
    switch,
    title,
)
import molecularnodes as mn

canvas = make_canvas()
protein = mn.Molecule.fetch("2K39")
align_frames(protein, "name CA and resid 1-70")
protein.add_style("cartoon").add_style("ball_and_stick")
cartoon, sticks = list(protein.styles)
canvas.look_at(protein, viewpoint=(70, 0, 37), margin=0.15)

heading = title(protein, "Every residue has a story.", position=(0.065, 0.07), size=64)
subtitle = title(
    protein,
    "Ubiquitin · PDB 2K39 · 116 NMR conformers",
    position=(0.067, 0.15),
    size=26,
    color="#a3b3c7",
)
footer = title(
    protein,
    "Interpolated NMR conformers · model order is not physical time",
    position=(0.065, 0.93),
    size=21,
    color="#74869e",
)
helix_note = callout(
    protein,
    "resid 23-34 and name CA",
    "α helix|Residues 23–34",
    offset=(-260, 40),
    size=38,
    color=GOLD,
)
tail_note = callout(
    protein,
    "resid 71-76 and name CA",
    "C-terminal tail|Residues 71–76",
    offset=(220, -60),
    size=38,
    color=CYAN,
)
# residue labels at the Cα of residues 8, 44 and 70, ProteinMotion's label_residues
residue_labels = []
for resid, offset in ((8, (-270, -90)), (44, (520, -100)), (70, (400, 150))):
    residue = protein.universe.select_atoms(f"resid {resid} and name CA").residues[0]
    residue_labels.append(
        callout(
            protein,
            f"resid {resid} and name CA",
            f"{residue.resname.title()} {resid}",
            offset=offset,
            size=34,
        )
    )

with canvas.timeline() as t:
    t.set_visible(sticks, False)
    show(t, heading)  # PM: Write(title, stroke_width=1.8), run_time=2.2
    t.wait(2.2)
    show(t, subtitle, footer)
    t.wait(0.6)
    show(t, helix_note, tail_note)
    # PM: helix.highlight(style="sphere"), tail.highlight(style="atoms") - the
    # regions are recoloured instead
    t.play(
        t.recolor(protein, GOLD, selection="resid 23-34"),
        t.recolor(protein, CYAN, selection="resid 71-76"),
        run_time=0.8,
    )
    t.wait(1.2)
    t.play(protein.play(0, 16), canvas.camera.orbit(20), run_time=4, easing="smooth")
    show(t, helix_note, tail_note, visible=False)
    t.wait(1.2)
    switch(t, [cartoon, sticks], [sticks])
    t.wait(0.3)
    show(t, *residue_labels)
    t.play(t.recolor(protein, GOLD, selection="resid 8 44 70"), run_time=0.8)
    t.wait(1.2)
    t.play(protein.play(16, 36), canvas.camera.orbit(23), run_time=5, easing="smooth")
    t.wait(1)

render(canvas, t, "labels_and_callouts")
