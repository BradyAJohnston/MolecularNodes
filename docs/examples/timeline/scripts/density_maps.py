"""1UBQ electron density from PDBe, with an animated contour and a moving slice.

After ProteinMotion's ``examples/density_maps.py``. The map is the PDBe
``1ubq.ccp4`` (a full unit cell), loaded through Molecular Nodes' density
entity and styled with ``Density Style ISO Surface``, whose ``threshold`` the
timeline keys. The slice is a thin box (an empty) read by a ``Select Cube``
node: a ``Delete Geometry`` after the surface keeps only what is inside it, and
the timeline moves the box through the cell.
"""

import urllib.request
import bpy
import numpy as np
from _common import GOLD, make_canvas, render, rgba, show, title
import molecularnodes as mn
from molecularnodes.blender import utils
from molecularnodes.entities.density import Grids
from molecularnodes.nodes import geometry as mg

MAP_URL = "https://www.ebi.ac.uk/pdbe/coordinates/files/1ubq.ccp4"


def add_slice(entity, style) -> bpy.types.Object:
    """
    Keep only the part of ``style``'s geometry inside a box empty, returned
    at its full size: scale and move it to cut a section.
    """
    points = utils.evaluated_points(entity.object)
    low, high = points.min(axis=0), points.max(axis=0)
    box = bpy.data.objects.new("slice", None)
    box.empty_display_type = "CUBE"
    box.location = tuple((low + high) / 2)
    # a cube empty spans -1 to 1, so its scale is half the box
    box.scale = tuple((high - low) / 2 * 1.05)
    bpy.context.scene.collection.objects.link(box)

    tree = entity.modifier_node_tree
    output = style.node.outputs["Geometry"]
    targets = [link.to_socket for link in output.links]
    for link in list(output.links):
        tree.links.remove(link)
    builder = entity.tree
    builder.activate_tree()
    try:
        inside = mg.SelectCube(object=box)
    finally:
        builder.deactivate_tree()
    delete = tree.nodes.new("GeometryNodeDeleteGeometry")
    tree.links.new(output, delete.inputs["Geometry"])
    tree.links.new(inside.node.outputs["Inverted"], delete.inputs["Selection"])
    for socket in targets:
        tree.links.new(delete.outputs["Geometry"], socket)
    return box


canvas = make_canvas()
map_path = mn.download.CACHE_DIR / "1ubq.ccp4"
if not map_path.exists():
    map_path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(MAP_URL, map_path)

protein = mn.Molecule.fetch("1UBQ")
# PM: protein.set_residue_opacity(0.08); helix.set_opacity(1). The Default
# material reads the Color attribute's alpha, so colour the context faint.
protein.set_color(rgba("#9fb3c8", alpha=0.08), "protein")
protein.set_color(rgba(GOLD), "resid 23-34")
protein.add_style("ball_and_stick", selection="protein")

density = Grids.load(map_path, style="density_iso_surface")
shell = density.styles[0]
# PM contours in sigma of the original map; the node's threshold is in map
# units, so convert from the grid's statistics
values = density.grid.grid
sigma = lambda n: float(values.mean() + n * values.std())  # noqa: E731
shell.i.threshold.socket.default_value = sigma(1.5)
box = add_slice(density, shell)
full = tuple(box.scale)
# PM: density.crop(helix, padding=2.5) - no crop yet, the whole cell is shown.

helix = protein.get_view("resid 23-34")
canvas.look_at(helix, viewpoint=(70, 0, 25), margin=0.55)
heading = title(protein, "Electron density", size=42)
detail = title(protein, "1UBQ · helix 23–34 · PDBe map", position=(0.06, 0.13), size=25)

z_low = box.location.z - full[2]
z_span = 2 * full[2]
with canvas.timeline() as t:
    t.set_visible(shell, False)
    show(t, heading, detail)
    # PM: FadeIn(shell)
    t.play(
        [
            t.fade(shell, 1.0),
            t.tween(shell.i.threshold, sigma(2.5)),
            canvas.camera.orbit(14),
        ],
        2,
    )
    t.play(t.tween(shell.i.threshold, sigma(1.5)), 1.5)
    # PM: local.slice("z", 0.15, ...) then section.animate.set_slice(0.85)
    centre = np.array(box.location)
    t.set((box, "scale"), (full[0], full[1], 0.03 * full[2]))
    t.set((box, "location"), (centre[0], centre[1], z_low + 0.15 * z_span))
    t.play(
        t.tween((box, "location"), (centre[0], centre[1], z_low + 0.85 * z_span)),
        3,
    )
    t.set((box, "scale"), full)
    t.wait(0.5)

render(canvas, t, "density_maps")
