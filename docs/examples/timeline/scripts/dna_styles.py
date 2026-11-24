"""DNA base styles, residue colours, transparency, labels and a surface.

After ProteinMotion's ``examples/dna_styles.py`` (1BNA).
"""

from _common import GOLD, callout, make_canvas, render, show, switch, title
import molecularnodes as mn

canvas = make_canvas()
dna = mn.Molecule.fetch("1BNA")
region = "nucleic and chainID A and resid 5-8"
other = "nucleic and chainID B"
# PM's base styles: slabs -> Rectangle, sticks -> Cylinder; "rings" and
# "ladder" have no cartoon equivalent yet (Style Ribbon draws base geometry)
dna.add_style("cartoon", selection=f"nucleic and not ({other})", base_shape="Rectangle")
dna.add_style("cartoon", selection=other, base_shape="Rectangle")
dna.add_style("cartoon", selection="nucleic", base_shape="Cylinder")
dna.add_style("ball_and_stick", selection="nucleic")
dna.add_style("surface", selection="nucleic")
slabs, slabs_b, sticks_bases, ball_stick, surface = list(dna.styles)
canvas.look_at(dna, viewpoint=(75, 0, -15), margin=0.05)
heading = title(dna, "DNA · 1BNA", position=(0.055, 0.88), size=30)
captions = {
    name: title(dna, name, position=(0.055, 0.10), size=46)
    for name in ("Base slabs", "Base sticks", "Residues 5–8", "Molecular surface")
}
note = callout(
    dna,
    region,
    "Residues 5–8|Color and opacity follow the selection",
    offset=(200, -60),
    size=30,
)
all_styles = [slabs, slabs_b, sticks_bases, ball_stick, surface]

with canvas.timeline() as t:
    t.set_visible([sticks_bases, ball_stick, surface], False)
    t.play(canvas.camera.focus(dna.get_view(region), fstop=5.6), 0)
    show(t, heading, captions["Base slabs"])
    t.play(canvas.camera.orbit(7), 1)
    t.play(canvas.camera.orbit(10), 2)
    show(t, captions["Base slabs"], visible=False)
    switch(t, all_styles, [sticks_bases])
    t.play(canvas.camera.orbit(6), 0.8)
    show(t, captions["Base sticks"])
    t.play(canvas.camera.orbit(7), 0.8)
    t.play(canvas.camera.orbit(10), 1.4)
    show(t, captions["Base sticks"], visible=False)
    switch(t, all_styles, [slabs, slabs_b])
    # PM: Colorize(region, gold, residue_delay=0.12), SetOpacity(other, 0.16).
    # The styles shown later are recoloured too, so the region stays gold.
    t.play(
        [
            t.recolor([slabs, ball_stick, surface], GOLD, selection=region),
            t.fade(slabs_b, 0.16),
        ],
        1.5,
    )
    show(t, note)
    t.play(canvas.camera.orbit(7), 1)
    switch(t, all_styles, [ball_stick])
    t.play(canvas.camera.orbit(12), 2)
    show(t, note, visible=False)
    switch(t, all_styles, [surface])
    show(t, captions["Molecular surface"])
    t.play(canvas.camera.orbit(17), 2.5)
    t.wait(0.5)

render(canvas, t, "dna_styles")
