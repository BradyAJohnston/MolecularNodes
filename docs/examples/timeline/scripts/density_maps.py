"""1UBQ electron density from PDBe, with an animated contour and a moving slice.

After ProteinMotion's ``examples/density_maps.py``. The map is the PDBe
``1ubq.ccp4`` (a full unit cell), loaded through Molecular Nodes' density
entity and styled with ``Density Style ISO Surface``, whose ``threshold``,
``slice_center`` and ``slice_width`` inputs the timeline keys.
"""

import urllib.request
from _common import GOLD, make_canvas, paint, render, show, title
import molecularnodes as mn
from molecularnodes.entities.density import Grids

MAP_URL = "https://www.ebi.ac.uk/pdbe/coordinates/files/1ubq.ccp4"

canvas = make_canvas()
map_path = mn.download.CACHE_DIR / "1ubq.ccp4"
if not map_path.exists():
    map_path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(MAP_URL, map_path)

protein = mn.Molecule.fetch("1UBQ")
# PM: protein.set_residue_opacity(0.08); helix.set_opacity(1). The Default
# material reads the Color attribute's alpha, so paint the context faint.
paint(protein, "protein", "#9fb3c8", alpha=0.08)
paint(protein, "resid 23-34", GOLD, alpha=1.0)
protein.add_style("ball_and_stick", selection="protein")

density = Grids.load(map_path, style="density_iso_surface")
shell = density.styles[0]
# PM contours in sigma of the original map; the node's threshold is in map
# units, so convert from the grid's statistics
values = density.grid.grid
sigma = lambda n: float(values.mean() + n * values.std())  # noqa: E731
shell.i.threshold.socket.default_value = sigma(1.5)
# PM: density.crop(helix, padding=2.5) - no crop yet, the whole cell is shown.

helix = protein.get_view("resid 23-34")
canvas.look_at(helix, viewpoint=(70, 0, 25), margin=0.55)
heading = title(protein, "Electron density", size=42)
detail = title(protein, "1UBQ · helix 23–34 · PDBe map", position=(0.06, 0.13), size=25)

with canvas.timeline() as t:
    show(t, heading, detail)
    # PM: FadeIn(shell) - no opacity channel on the density style yet
    t.play(t.tween(shell.i.threshold, sigma(2.5)), canvas.camera.orbit(14), run_time=2)
    t.play(t.tween(shell.i.threshold, sigma(1.5)), run_time=1.5)
    # PM: local.slice("z", 0.15, ...) then section.animate.set_slice(0.85): the
    # ISO surface node has a slice built in
    # ISO surface node; width and centre are fractions of the box per axis
    t.set(shell.i.slice_width, (0.5, 0.5, 0.03))
    t.set(shell.i.slice_center, (0.5, 0.5, 0.15))
    t.play(t.tween(shell.i.slice_center, (0.5, 0.5, 0.85)), run_time=3)
    t.set(shell.i.slice_width, (0.5, 0.5, 0.5))
    t.wait(0.5)

render(canvas, t, "density_maps")
