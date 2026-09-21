"""A continuous calmodulin feature tour, with no scene cuts.

After ProteinMotion's ``examples/feature_showcase.py``. Chapters that need a
node Molecular Nodes does not have yet (the contact-guided morph into troponin
C, hydrogen bonds, electrostatics, write-on text) are marked in comments; the
NMR playback uses the deposited 1CFC ensemble directly rather than a mapped
XTC.
"""

from _common import (
    BLUE,
    GOLD,
    TEAL,
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

MUTED = "#9eafc4"

canvas = make_canvas()
p = mn.Molecule.fetch("1CLL")
helix, lobe = "protein and resid 69-88", "protein and resid 5-38"
paint(p, None, "#8299ad")
context = context_material(0.0)
p.add_style("cartoon", selection=f"{helix} or {lobe}")
p.add_style(
    "cartoon", selection=f"protein and not ({helix} or {lobe})", material=context
)
p.add_style("ribbon", selection="protein")
p.add_style("ball_and_stick", selection="protein")
p.add_style("surface", selection="protein")
cartoon, cartoon_ctx, ribbon, sticks, surface = list(p.styles)
styles = [cartoon, cartoon_ctx, ribbon, sticks, surface]
hide(ribbon, sticks, surface)
canvas.look_at(p, viewpoint=(72, 0, 12), margin=0.35)

chapters = {}
for n, (name, head, detail) in enumerate(
    [
        (
            "Representations",
            "One protein. Every view.",
            "Cartoon · ribbon · ball and stick · surface",
        ),
        (
            "Color + transparency",
            "Color flows through residues.",
            "One selection across every representation",
        ),
        (
            "Focus + annotation",
            "Stay with the same structure.",
            "Eased focus · callouts",
        ),
        (
            "States",
            "Calmodulin keeps moving.",
            "1CFC NMR conformers · illustrative interpolation, not MD",
        ),
        (
            "One continuous scene",
            "Every detail. One connected story.",
            "Authored in Python · rendered in Blender",
        ),
    ],
    start=1,
):
    chapters[name] = (
        title(
            p, f"{n:02d} / {name.upper()}", position=(0.06, 0.055), size=24, color=TEAL
        ),
        title(p, head, position=(0.06, 0.103), size=62),
        title(p, detail, position=(0.06, 0.925), size=23, color=MUTED),
    )
identity = title(p, "CALMODULIN · 1CLL", position=(0.72, 0.055), size=23, color=MUTED)
badges = {
    name: title(p, name, position=(0.72, 0.83), size=34, color=color)
    for name, color in (
        ("Cartoon", TEAL),
        ("Ribbon", BLUE),
        ("Ball and stick", GOLD),
        ("Molecular surface", TEAL),
    )
}
note = callout(
    p,
    helix,
    "Central helix|Calmodulin · residues 69–88",
    offset=(-260, 40),
    size=38,
    color=TEAL,
)
current = None


def heading(t, name):
    global current
    if current:
        show(t, *chapters[current], visible=False)
    show(t, *chapters[name])
    current = name


with canvas.timeline() as t:
    heading(t, "Representations")
    show(t, identity, badges["Cartoon"])
    t.play(canvas.camera.orbit(7), run_time=1.8)
    last = badges["Cartoon"]
    for style, name in (
        (ribbon, "Ribbon"),
        (sticks, "Ball and stick"),
        (surface, "Molecular surface"),
    ):
        show(t, last, visible=False)
        show_only(t, styles, [style])
        show(t, badges[name])
        t.play(canvas.camera.orbit(7), run_time=1.9)
        t.play(canvas.camera.orbit(5), run_time=1.1)
        last = badges[name]
    show(t, last, visible=False)
    show_only(t, styles, [cartoon, cartoon_ctx])

    heading(t, "Color + transparency")
    # PM: Colorize(helix, TEAL, residue_delay=0.065), Colorize(lobe, GOLD, ...):
    # staggered colour needs the planned Animate Stagger node; painted below
    paint(p, helix, TEAL)
    paint(p, lobe, GOLD)
    t.play(canvas.camera.orbit(7), run_time=2.5)
    t.play(
        t.tween(context.node.i.transparency, 0.8), canvas.camera.orbit(6), run_time=2
    )
    t.play(
        t.tween(context.node.i.transparency, 0.0), canvas.camera.orbit(5), run_time=1.5
    )

    heading(t, "Focus + annotation")
    # PM: helix.highlight(style="sphere"), highlight(style="box") - no highlight style yet
    t.play(
        canvas.camera.look_at(p.get_view(helix), margin=0.3),
        t.tween(context.node.i.transparency, 0.95),
        run_time=2.8,
    )
    show(t, note)
    t.play(canvas.camera.orbit(13), run_time=2.6)
    show(t, note, visible=False)
    t.play(
        canvas.camera.look_at(p, margin=0.35),
        t.tween(context.node.i.transparency, 0.0),
        run_time=2.5,
    )

    heading(t, "States")
    # PM plays a 1CFC ensemble mapped onto 1CLL's atoms, then Deform + Morph.
    # Here the deposited ensemble is a second entity, cut in at this point.
    t.wait(0.5)

render(canvas, t, "feature_showcase_a")

nmr = mn.Molecule.fetch("1CFC")
align_frames(nmr, "name CA")
paint(nmr, None, "#8299ad")
nmr.add_style("cartoon", selection="protein")
p.object.hide_render = True
canvas.look_at(nmr, viewpoint=(72, 0, 12), margin=0.35)
nmr_identity = title(
    nmr, "CALMODULIN · 1CFC", position=(0.72, 0.055), size=23, color=MUTED
)
with canvas.timeline(start=t.frame_end + 1) as t2:
    show(t2, nmr_identity)
    t2.play(
        nmr.play(0, nmr.universe.trajectory.n_frames - 1),
        canvas.camera.orbit(14),
        run_time=7,
    )
    # PM: BackboneMorph(p, troponin, match=...) - planned Morph To Attribute node;
    # HydrogenBonds, Distance, Electrostatics - planned analysis nodes
    t2.play(canvas.camera.orbit(10), run_time=2)
    t2.wait(1)

render(canvas, t2, "feature_showcase_b")
