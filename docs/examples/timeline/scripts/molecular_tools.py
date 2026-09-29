"""Residue styling, a surface, a live ruler, and the pieces that are not here yet.

After ProteinMotion's ``examples/molecular_tools.py`` (2K39). Hydrogen-bond
and electrostatic contact rulers need the planned ``Hydrogen Bonds`` and
``Coulomb Contacts`` nodes; the Cα distance ruler is the ``com_distance``
annotation, which follows the ensemble as it plays.
"""

from _common import (
    CYAN,
    GOLD,
    INK,
    align_frames,
    callout,
    context_material,
    hide,
    make_canvas,
    paint,
    render,
    show,
    show_only,
    title,
)
import molecularnodes as mn

canvas = make_canvas()
protein = mn.Molecule.fetch("2K39")
align_frames(protein, "name CA and resid 1-70")
helix, strand = "resid 23-34", "resid 2-7"
rest = f"protein and not ({helix} or {strand})"
paint(protein, None, "#8299ad")
paint(protein, helix, CYAN)
paint(protein, strand, GOLD)
context = context_material(0.0)
protein.add_style("cartoon", selection=f"{helix} or {strand}")
protein.add_style("cartoon", selection=rest, material=context)
protein.add_style("ball_and_stick", selection="protein")
protein.add_style("surface", selection="protein")
cartoon, faint, sticks, surface = list(protein.styles)
hide(sticks, surface)
canvas.look_at(protein, viewpoint=(70, 0, 35), margin=0.25)

heading = title(protein, "Color is part of the story.", size=64)
note = title(
    protein,
    "Residue colors · local transparency",
    position=(0.065, 0.155),
    size=26,
    color=INK,
)
footer = title(
    protein,
    "Ubiquitin · 2K39 · interpolated conformers are not physical time",
    position=(0.065, 0.93),
    size=21,
    color="#74869e",
)
helix_note = callout(
    protein, helix, "α helix|Residues 23–34", offset=(-260, 40), size=36, color=CYAN
)
strand_note = callout(
    protein, strand, "β strand|Residues 2–7", offset=(220, -60), size=36, color=GOLD
)
ruler = protein.annotations.add_com_distance(
    selection1="resid 23 and name CA",
    selection2="resid 34 and name CA",
    text1="Cα 23",
    text2="Cα 34",
)
ruler.visible = False
surface_note = title(
    protein,
    "Surface · meshes follow the coordinates",
    position=(0.065, 0.155),
    size=26,
    color=INK,
)

with canvas.timeline() as t:
    show(t, heading, note, footer)
    # PM: Colorize(helix, ..., residue_delay=0.09) - staggered colour needs the
    # planned Animate Stagger node; the regions are painted from the start
    t.wait(2.5)
    show(t, helix_note, strand_note)
    t.wait(1.2)
    t.play(
        t.tween(context.node.i.transparency, 0.88), run_time=1.3
    )  # SetOpacity(rest, 0.12)
    t.play(canvas.camera.orbit(20), run_time=2)
    t.play(t.tween(context.node.i.transparency, 0.0), run_time=1)
    show(t, helix_note, strand_note, visible=False)
    show_only(t, [cartoon, faint, sticks, surface], [sticks])
    show(t, ruler)
    t.wait(1.5)
    t.play(protein.play(0, 12), run_time=3, easing="smooth")
    show(t, ruler, note, visible=False)
    show_only(t, [cartoon, faint, sticks, surface], [surface])
    show(t, surface_note)
    t.play(protein.play(12, 24), canvas.camera.orbit(23), run_time=3, easing="smooth")
    t.play(canvas.camera.orbit(23), run_time=2)
    t.wait(0.6)

render(canvas, t, "molecular_tools_styling")

# PM's second scene, InteractionsAndDistances: hydrogen-bond rulers with live
# donor-acceptor distances and screened electrostatic contacts. Both need
# analysis nodes that do not exist yet:
#
#   hbonds = protein.measure.hbonds(donors="resid 23-34 and name N", max_distance=3.5, min_angle=150)
#   t.play(hbonds.show(mode="3d", show_distances=True), run_time=1.6)
#   field = protein.measure.contacts(charges="formal", dielectric=80, screening_length=8)
#   t.play(field.show(mode="2d", max_pairs=3), run_time=1.6)
