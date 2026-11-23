"""Calmodulin in focus: a continuous camera study of PDB 1CLL, chain A.

After ProteinMotion's ``examples/calmodulin_in_focus.py``, a 68-second film.
Coordinates stay fixed; the camera, the depth of field and the cross-fades
between representations are the story. The backbone hydrogen-bond network and the colour legend
are not here yet.
"""

from _common import (
    BLUE,
    GOLD,
    INK,
    TEAL,
    WHITE,
    callout,
    make_canvas,
    render,
    show,
    title,
)
import molecularnodes as mn
from molecularnodes.nodes import geometry as mg

canvas = make_canvas()
p = mn.Molecule.fetch("1CLL")
helix, far_helix = "protein and resid 5-19", "protein and resid 118-128"
backbone = f"{helix} and backbone"
p.set_color(BLUE)
p.set_color(GOLD, helix)
# the region and its context are separate styles, so the context can be faded
p.add_style("cartoon", selection=f"{helix} or {far_helix}")
p.add_style("cartoon", selection=f"protein and not ({helix} or {far_helix})")
p.add_style("ball_and_stick", selection=backbone)
p.add_style("ball_and_stick", selection=f"protein and not ({backbone})")
p.add_style(
    "surface",
    selection="protein",
    color=lambda: mg.ColorAttributeMap(name="b_factor", min=0, max=50),
)
p.add_style("ribbon", selection="protein")
cartoon, cartoon_ctx, sticks, sticks_ctx, surface, ribbon = list(p.styles)
canvas.look_at(p, viewpoint=(72, 0, -30), margin=0.18)

captions = {}
for key, (head, detail, color) in {
    "open": ("Calmodulin", "1CLL · a study in depth and focus", WHITE),
    "helix": ("An alpha helix", "Chain A · residues 5–19", GOLD),
    "inside": ("Inside the helix", "Backbone N···O hydrogen-bond candidates", WHITE),
    "moving": ("Moving the focus", "Gold: 5–19   /   Teal: 118–128", WHITE),
    "focus_a": (
        "Focus · residues 5–19",
        "The lens follows the selected Cα atoms",
        GOLD,
    ),
    "focus_b": (
        "Focus · residues 118–128",
        "A continuous focus pull to the second helix",
        TEAL,
    ),
    "surface": (
        "The molecular surface",
        "Color follows the deposited B factors",
        WHITE,
    ),
    "back": ("Back to the backbone", "One structure · continuous motion", WHITE),
    "end": ("Calmodulin in focus", "Made with Molecular Nodes and Blender", WHITE),
}.items():
    captions[key] = (
        title(p, head, position=(0.06, 0.07), size=48, color=color),
        title(p, detail, position=(0.06, 0.135), size=25, color=INK),
    )
identity = title(
    p, "Molecular Nodes / EEVEE", position=(0.78, 0.06), size=22, color=INK
)
note = callout(p, helix, "Residues 5–19", offset=(220, 60), size=30)
current = None


def caption(t, key):
    global current
    if current is not None:
        show(t, *captions[current], visible=False)
    show(t, *captions[key])
    current = key


ca = lambda sel: p.get_view(f"{sel} and name CA")  # noqa: E731

with canvas.timeline() as t:
    t.set_visible([sticks, sticks_ctx, surface, ribbon], False)
    # 00-08: a wide reveal and an unhurried orbit
    t.play(canvas.camera.focus(ca(helix), fstop=5.6), run_time=0)
    caption(t, "open")
    show(t, identity)
    t.wait(2)
    t.play(canvas.camera.orbit(17), canvas.camera.zoom(54), run_time=6)
    # 08-17: isolate one helix, keeping a translucent context
    caption(t, "helix")
    t.play(
        canvas.camera.look_at(p.get_view(helix), margin=0.3),
        canvas.camera.focus(ca(helix), fstop=3.5),
        t.fade(cartoon_ctx, 0.1),
        run_time=3,
    )
    show(t, note)
    t.play(canvas.camera.orbit(5), run_time=1.5)
    t.play(canvas.camera.orbit(16), run_time=4.5)
    # 17-27: into an atomic view of the backbone
    caption(t, "inside")
    show(t, note, visible=False)
    t.play(
        t.fade([cartoon, cartoon_ctx], 0.0),
        t.fade(sticks, 1.0),
        t.fade(sticks_ctx, 0.035),
        canvas.camera.focus(ca(backbone), fstop=12),
        canvas.camera.look_at(p.get_view(backbone), margin=0.05),
        run_time=3,
    )
    # PM: HydrogenBonds(...).highlight(mode="3d") - see alpha_helix_hbonds.py
    t.play(canvas.camera.orbit(6), run_time=2)
    t.play(
        canvas.camera.orbit(24),
        canvas.camera.focus(ca("resid 11-14"), fstop=10),
        run_time=5,
    )
    # 27-35: back out and establish two focus targets
    caption(t, "moving")
    t.play(
        t.fade([sticks, sticks_ctx], 0.0),
        t.fade([cartoon, cartoon_ctx], 1.0),
        canvas.camera.look_at(p.get_view("protein"), margin=0.12),
        canvas.camera.focus(ca(helix), fstop=3.2),
        run_time=4,
    )
    t.play(t.fade(cartoon_ctx, 0.24), canvas.camera.orbit(17), run_time=4)
    # 35-45: hold the framing so the focus transfer is visible
    caption(t, "focus_a")
    t.play(
        canvas.camera.focus(ca(helix), fstop=1.4), canvas.camera.orbit(6), run_time=5
    )
    caption(t, "focus_b")
    t.play(
        canvas.camera.focus(ca(far_helix), fstop=1.4),
        canvas.camera.orbit(6),
        run_time=5,
    )
    # 45-55: the surface, coloured by B factor
    caption(t, "surface")
    t.play(
        t.fade([cartoon, cartoon_ctx], 0.0),
        t.fade(surface, 1.0),
        canvas.camera.look_at(p.get_view("protein"), margin=0.15),
        canvas.camera.focus(ca(far_helix), fstop=5.6),
        run_time=3,
    )
    # PM: ColorLegend(scale, title="Cα B factor") - no legend annotation yet
    t.play(canvas.camera.orbit(32), run_time=7)
    # 55-68: back to the backbone, and a last focus pull
    caption(t, "back")
    t.play(
        t.fade(surface, 0.0),
        t.fade(ribbon, 1.0),
        canvas.camera.look_at(p.get_view("protein"), margin=0.16),
        canvas.camera.focus(ca(helix), fstop=3.5),
        run_time=4,
    )
    caption(t, "end")
    t.play(
        t.fade(ribbon, 0.0),
        t.fade([cartoon, cartoon_ctx], 1.0),
        canvas.camera.orbit(24),
        canvas.camera.zoom(54),
        canvas.camera.focus(ca(helix), fstop=2),
        run_time=6,
    )
    # PM: FadeOut(p)
    t.play(t.fade(p, 0.0), canvas.camera.orbit(3), run_time=3)

render(canvas, t, "calmodulin_in_focus", frames=8)
