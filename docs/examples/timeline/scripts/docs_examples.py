"""Small scenes, one per documentation snippet.

After ProteinMotion's ``examples/docs_examples.py``. Run all of them, or one:

    uv run python docs/examples/timeline/scripts/docs_examples.py out_dir [name]

Scenes whose subject Molecular Nodes cannot animate yet (write-on text, 3D
highlights, charge contacts, deformation) render what can be shown and say
what is missing in a comment.
"""

import sys
from _common import (
    CYAN,
    GOLD,
    align_frames,
    callout,
    design_pixels,
    make_canvas,
    render,
    rgba,
    show,
    spin,
    switch,
    title,
)
import molecularnodes as mn
from molecularnodes.nodes import geometry as mg

SCENES = {}


def scene(fn):
    SCENES[fn.__name__] = fn
    return fn


def ubiquitin(canvas, style="cartoon", viewpoint=(70, 0, 37), margin=0.05, **kwargs):
    p = mn.Molecule.fetch("1UBQ")
    p.add_style(style, selection="protein", **kwargs)
    canvas.look_at(p, viewpoint=viewpoint, margin=margin)
    return p


@scene
def motion(canvas):
    p = ubiquitin(canvas)
    p.add_style("ribbon", selection="protein")
    cartoon, ribbon = list(p.styles)
    with canvas.timeline() as t:
        t.set_visible(p, False)
        t.play(t.fade(cartoon, 1.0), run_time=0.8)  # PM: FadeIn(p)
        t.play(spin(p, 180), run_time=3)
        switch(t, [cartoon, ribbon], [ribbon])
        t.play(
            t.tween((p.object, "location"), (0.5, 0, 0)),
            canvas.camera.orbit(34),
            run_time=2,
        )
        t.wait(0.5)
    return t


@scene
def representations(canvas):
    p = ubiquitin(canvas)
    for name in ("ribbon", "ball_and_stick", "surface"):
        p.add_style(name, selection="protein")
    styles = list(p.styles)
    with canvas.timeline() as t:
        t.set_visible(styles[1:], False)
        t.wait(1)
        for style in styles[1:]:
            switch(t, styles, [style])
            t.play(spin(p, 26), run_time=1.5)
    return t


@scene
def region_focus(canvas):
    p = ubiquitin(canvas)
    p.set_color(GOLD, "resid 23-34")  # PM: helix.highlight(style="box")
    helix = p.get_view("resid 23-34")
    with canvas.timeline() as t:
        t.play(canvas.camera.look_at(helix, margin=0.3), run_time=1.5)
        t.play(canvas.camera.orbit(23), run_time=2)
        t.play(canvas.camera.look_at(p, margin=0.05), run_time=1.5)
    return t


@scene
def highlights(canvas):
    # PM: helix.highlight(style="sphere" | "box" | "atoms") - the planned Style
    # Highlight node. The helix is coloured and shown with each representation.
    p = ubiquitin(canvas)
    p.set_color(GOLD, "resid 23-34")
    with canvas.timeline() as t:
        t.play(spin(p, 40), run_time=3)
    return t


@scene
def writing(canvas):
    # PM: Write(title, lag_ratio=0.12) / Unwrite(title) - vector write-on needs a
    # ``progress`` property on text annotations; the title is cut in and out.
    p = ubiquitin(canvas)
    heading = title(p, "Protein motion", position=(0.1, 0.38), size=130)
    with canvas.timeline() as t:
        show(t, heading)
        t.wait(2.5)
        t.wait(1)
        show(t, heading, visible=False)
        t.wait(0.5)
    return t


@scene
def text_placement(canvas):
    p = ubiquitin(canvas)
    heading = title(p, "α helix / β sheet", position=(0.06, 0.08), size=80, color=CYAN)
    footer = title(p, "PDB 1UBQ", position=(0.8, 0.9), size=40)
    with canvas.timeline() as t:
        show(t, heading, footer)
        t.wait(1)
        t.play(t.tween((heading, "location"), (0.10, 0.65)), run_time=1)
        # PM: title.animate.set_opacity(0.4) - labels have no opacity; shrink instead
        t.play(t.tween((heading, "text_size"), design_pixels(50)), run_time=0.5)
        t.wait(1)
    return t


@scene
def callout_scene(canvas):
    p = ubiquitin(canvas, margin=0.05)
    note = callout(
        p,
        "resid 23-34 and name CA",
        "α helix|Residues 23–34",
        offset=(-260, 40),
        size=48,
    )
    p.set_color(GOLD, "resid 23-34")  # PM: helix.highlight(style="box")
    with canvas.timeline() as t:
        show(t, note)
        t.wait(1)
        t.play(canvas.camera.orbit(29), run_time=3)
        t.wait(1)
    return t


@scene
def residue_labels(canvas):
    p = ubiquitin(canvas, margin=0.05)
    labels = []
    for resid, offset in ((8, (-270, -90)), (44, (440, -100)), (70, (380, 150))):
        residue = p.universe.select_atoms(f"resid {resid} and name CA").residues[0]
        labels.append(
            callout(
                p,
                f"resid {resid} and name CA",
                f"{residue.resname.title()} {resid}",
                offset=offset,
                size=40,
            )
        )
    with canvas.timeline() as t:
        show(t, *labels)
        t.wait(1)
        t.play(canvas.camera.orbit(20), run_time=2)
        t.wait(1)
    return t


@scene
def residue_colors(canvas):
    p = mn.Molecule.fetch("1UBQ")
    p.set_color(CYAN, "resid 23-34")
    p.set_color(GOLD, "resid 2-7")
    # opacity through the Color alpha, which the Default material reads
    p.set_color(rgba("#9fb3c8", alpha=0.2), "resid 71-76")
    p.add_style("cartoon", selection="protein")
    canvas.look_at(p, viewpoint=(70, 0, 37), margin=0.05)
    with canvas.timeline() as t:
        t.play(spin(p, 40), run_time=3)
    return t


@scene
def color_change(canvas):
    # PM: Colorize(helix, cyan, residue_delay=0.09); the colour blends in
    # without the per-residue stagger
    p = mn.Molecule.fetch("1UBQ")
    p.add_style("cartoon", selection="protein")
    canvas.look_at(p, viewpoint=(70, 0, 37), margin=0.05)
    with canvas.timeline() as t:
        t.wait(1)
        t.play(
            t.recolor(p, CYAN, selection="resid 23-34"),
            t.recolor(p, GOLD, selection="resid 2-7"),
            run_time=1.5,
        )
        t.wait(1)
    return t


@scene
def opacity(canvas):
    p = mn.Molecule.fetch("1UBQ")
    p.set_color(CYAN, "resid 23-34")
    p.add_style("ball_and_stick", selection="resid 23-34")
    p.add_style("ball_and_stick", selection="protein and not resid 23-34")
    helix = p.styles[0]
    canvas.look_at(p, viewpoint=(70, 0, 37), margin=0.05)
    with canvas.timeline() as t:
        t.wait(0.5)
        t.play(t.fade(helix, 0.1), run_time=2)  # SetOpacity(helix, 0.1)
        t.wait(0.7)
        t.play(t.fade(helix, 1.0), run_time=1.5)
        t.wait(0.5)
    return t


@scene
def surface(canvas):
    p = ubiquitin(canvas, style="surface")
    with canvas.timeline() as t:
        t.play(spin(p, 46), run_time=3)
    return t


@scene
def distances(canvas):
    p = ubiquitin(canvas, style="ball_and_stick", margin=0.1)
    ruler = p.annotations.add_com_distance(
        selection1="resid 23 and name CA",
        selection2="resid 34 and name CA",
        text1="Cα 23",
        text2="Cα 34",
    )
    ruler.visible = False
    heading = title(p, "Distance", size=48)
    with canvas.timeline() as t:
        show(t, ruler, heading)
        t.wait(1.5)
        t.play(canvas.camera.orbit(20), run_time=2)
        show(t, ruler, heading, visible=False)
        t.wait(0.5)
    return t


@scene
def hydrogen_bonds(canvas):
    # PM: HydrogenBonds(p, donors=..., max_distance=3.5, min_angle=150).highlight()
    # MDAnalysis needs explicit hydrogens, so this is the 2K39 NMR model of
    # ubiquitin rather than the 1UBQ crystal structure
    p = mn.Molecule.fetch("2K39")
    p.add_style("cartoon", selection="protein")
    hbonds = p.calculations.hbonds(
        selection="protein", d_a_cutoff=3.5, d_h_a_angle_cutoff=150
    )
    with p.tree as tree:
        (
            tree.atoms
            >> mg.StyleHBonds(
                data_object=hbonds.object,
                radius=0.2,
                color=rgba(GOLD),
                material=mn.material.append_material("Default"),
            )
            >> tree.join
        )
    network = p.styles["Style HBonds"]
    canvas.look_at(p, viewpoint=(70, 0, 37), margin=0.1)
    with canvas.timeline() as t:
        t.set_visible(network, False)
        t.wait(0.5)
        t.play(t.fade(network, 1.0), run_time=1)
        t.wait(0.5)
        t.play(canvas.camera.orbit(34), run_time=3)
        t.wait(1)
    return t


@scene
def charge_contacts(canvas):
    # PM: Electrostatics(p, charges="formal", ...).highlight(mode="2d") needs the
    # planned Coulomb Contacts node. Nothing to animate yet; a still turntable.
    ubiquitin(canvas, style="ball_and_stick")
    with canvas.timeline() as t:
        t.play(canvas.camera.orbit(29), run_time=3)
        t.wait(1)
    return t


@scene
def deformation(canvas):
    # PM: Deform(p, stretch) then Morph(p, rest) needs the planned Morph To
    # Attribute node. Scaling the object stands in for the stretch.
    p = ubiquitin(canvas, margin=0.05)
    with canvas.timeline() as t:
        t.play(t.tween((p.object, "scale"), (1.4, 1, 1)), run_time=2)
        t.wait(0.5)
        t.play(t.tween((p.object, "scale"), (1, 1, 1)), run_time=2)
        t.wait(0.5)
    return t


@scene
def nmr_states(canvas):
    p = mn.Molecule.fetch("2K39")
    align_frames(p, "name CA and resid 1-70")
    p.add_style("cartoon", selection="protein")
    canvas.look_at(p, viewpoint=(70, 0, 37), margin=0.05)
    with canvas.timeline() as t:
        t.play(p.play(0, 2), run_time=6, easing="smooth")
        t.wait(0.5)
    return t


if __name__ == "__main__":
    names = sys.argv[2:] or list(SCENES)
    for name in names:
        canvas = make_canvas()
        canvas.clear()
        t = SCENES[name](canvas)
        render(canvas, t, f"docs_{name}")
