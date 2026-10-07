"""A ubiquitin NMR ensemble (2K39), aligned on the core and coloured by RMSF.

After ProteinMotion's ``examples/synchronized_plots.py``. The contact map,
sequence track, distance trace and legend are plot annotations Molecular
Nodes does not have yet (comparison doc, 5.4); the playback and the RMSF
colouring are here.
"""

import numpy as np
from _common import align_frames, make_canvas, render, show, title
from MDAnalysis.analysis import rms
import molecularnodes as mn
from molecularnodes.nodes import geometry as mg

canvas = make_canvas()
protein = mn.Molecule.fetch("2K39")
align_frames(protein, "name CA and resid 1-70")
# per-residue RMSF over the aligned ensemble, stored as an attribute
u = protein.universe
ca = u.select_atoms("name CA")
rmsf = rms.RMSF(ca).run().results.rmsf
per_atom = np.zeros(len(u.atoms), dtype=np.float32)
for atom, value in zip(ca, rmsf):
    per_atom[atom.residue.atoms.ix] = value
protein.store_named_attribute(per_atom, "rmsf", atype="FLOAT")
protein.add_style(
    "cartoon", color=lambda: mg.ColorAttributeMap(name="rmsf", min=0, max=8)
)
# PM: selected.highlight(style="box") on residues 23-34 - no highlight style yet.
canvas.look_at(protein, viewpoint=(70, 0, 40), margin=0.2)
heading = title(
    protein, "An NMR ensemble with synchronized plots", position=(0.05, 0.05), size=35
)
detail = title(protein, "2K39 · deposited model order", position=(0.05, 0.105), size=24)
# PM: ContactMap(...), TimeSeriesPlot.distance(...), SequenceTrack(...), ColorLegend(...)

with canvas.timeline() as t:
    show(t, heading, detail)
    t.wait(0.5)
    t.play(protein.play(0, u.trajectory.n_frames - 1), run_time=9, easing="linear")
    t.wait(0.5)

render(canvas, t, "synchronized_plots")
