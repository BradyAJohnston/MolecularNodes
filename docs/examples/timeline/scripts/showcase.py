"""Ubiquitin: representations, a camera sweep and a synthetic trajectory.

After ProteinMotion's ``examples/showcase.py``. The animated motion is
synthetic (a sine-wave displacement), not an MD simulation; it is built as an
in-memory MDAnalysis trajectory so the entity plays it like any other.
"""

import MDAnalysis as mda
import numpy as np
from _common import hide, make_canvas, render, show_only, spin, world_points
from MDAnalysis.coordinates.memory import MemoryReader
import molecularnodes as mn


def synthetic_motion(xyz, phase):
    centered = xyz - xyz.mean(0)
    x, y, z = centered.T
    displacement = np.column_stack(
        (
            1.6 * np.sin(y * 0.12 + phase),
            0.9 * np.sin(z * 0.13 + phase),
            1.6 * np.sin(x * 0.10 + phase),
        )
    )
    return xyz + displacement


canvas = make_canvas()
pdb = mn.download.StructureDownloader().download("1UBQ", format="pdb")
u = mda.Universe(pdb)
original = u.atoms.positions.copy()
frames = [synthetic_motion(original, phase) for phase in np.linspace(0, 2 * np.pi, 61)]
u.load_new(np.stack(frames), format=MemoryReader)
protein = mn.Molecule(u)
protein.add_style("cartoon", selection="protein").add_style(
    "ribbon", selection="protein"
).add_style("ball_and_stick", selection="protein")
cartoon, ribbon, sticks = list(protein.styles)
hide(ribbon, sticks)
canvas.look_at(protein, viewpoint=(75, 0, -20), margin=0.05)

with canvas.timeline() as t:
    # PM: FadeIn(p) - no opacity channel yet
    t.play(spin(protein, 126), run_time=3.2)
    show_only(t, [cartoon, ribbon, sticks], [ribbon])
    t.play(spin(protein, 50, axis="x"), run_time=2)
    show_only(t, [cartoon, ribbon, sticks], [sticks])
    t.play(canvas.camera.orbit(50), run_time=2)
    show_only(t, [cartoon, ribbon, sticks], [cartoon])
    # PM: Deform(p, stretch) then Morph(p, original) - needs Morph To Attribute
    t.play(protein.play(0, 60), spin(protein, 85), run_time=4, easing="linear")
    t.wait(0.7)

render(canvas, t, "showcase")

# Gallery: three copies side by side, each spinning
canvas.clear()
copies = []
for offset, style in ((-3.6, "cartoon"), (0.0, "ribbon"), (3.6, "ball_and_stick")):
    mol = mn.Molecule.fetch("1UBQ")
    mol.add_style(style)
    mol.object.location = (offset, 0.0, 0.0)
    copies.append(mol)
canvas.look_at(world_points(*copies), viewpoint="front", margin=0.1)
with canvas.timeline() as t:
    t.play(*(spin(mol, 360) for mol in copies), run_time=12, easing="linear")

render(canvas, t, "showcase_gallery")
