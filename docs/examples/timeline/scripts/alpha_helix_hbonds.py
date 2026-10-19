"""A controlled alpha-helix hydrogen-bond test built from ideal geometry.

After ProteinMotion's ``examples/alpha_helix_hbonds.py``. The helix is built
in numpy (phi -57, psi -47, omega 180, explicit amide H) into an MDAnalysis
Universe. The i -> i+4 O···N pairs are drawn with ``com_distance`` annotations;
a geometric hydrogen-bond detector (the planned ``Hydrogen Bonds`` node) is
what ProteinMotion animates here.
"""

import MDAnalysis as mda
import numpy as np
from _common import GOLD, context_material, make_canvas, render, show, title
import molecularnodes as mn


def unit(v):
    return v / np.linalg.norm(v)


def place(a, b, c, length, angle, torsion):
    u = unit(c - b)
    normal = unit(np.cross(b - a, u))
    side = np.cross(normal, u)
    theta, tau = np.deg2rad([angle, torsion])
    return c + length * (
        -np.cos(theta) * u + np.sin(theta) * (np.cos(tau) * side + np.sin(tau) * normal)
    )


def ideal_helix(n=16, phi=-57, psi=-47) -> mda.Universe:
    bb = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.458, 0.0, 0.0]),
        np.array(
            [
                1.458 + 1.525 * np.cos(np.deg2rad(68.8)),
                1.525 * np.sin(np.deg2rad(68.8)),
                0.0,
            ]
        ),
    ]
    for _ in range(n):
        bb.append(place(*bb[-3:], 1.329, 116.2, psi))
        bb.append(place(*bb[-3:], 1.458, 121.7, 180))
        bb.append(place(*bb[-3:], 1.525, 111.2, phi))
    names, elements, resindex, xyz, bonds = [], [], [], [], []
    for j in range(n):
        nn, ca, c = bb[j * 3 : j * 3 + 3]
        nn1 = bb[j * 3 + 3]
        o = c - unit(unit(ca - c) + unit(nn1 - c)) * 1.231
        atoms = [("N", "N", nn), ("CA", "C", ca), ("C", "C", c), ("O", "O", o)]
        if j:
            h = nn - unit(unit(bb[j * 3 - 1] - nn) + unit(ca - nn)) * 1.01
            atoms.append(("H", "H", h))
        start = len(names)
        for name, element, position in atoms:
            names.append(name)
            elements.append(element)
            resindex.append(j)
            xyz.append(position)
        bonds += [(start, start + 1), (start + 1, start + 2), (start + 2, start + 3)]
        if j:
            bonds += [
                (start - 4 + 2 if j == 1 else start - 5 + 2, start),
                (start, start + 4),
            ]
    u = mda.Universe.empty(
        len(names),
        n_residues=n,
        atom_resindex=resindex,
        residue_segindex=[0] * n,
        trajectory=True,
    )
    u.add_TopologyAttr("name", names)
    u.add_TopologyAttr("type", elements)
    u.add_TopologyAttr("element", elements)
    # the distance annotations use centres of mass, so masses are needed
    masses = {"N": 14.007, "C": 12.011, "O": 15.999, "H": 1.008}
    u.add_TopologyAttr("mass", [masses[e] for e in elements])
    u.add_TopologyAttr("resname", ["ALA"] * n)
    u.add_TopologyAttr("resid", list(range(1, n + 1)))
    u.add_TopologyAttr("chainID", ["A"] * len(names))
    u.add_TopologyAttr("segid", ["A"])
    u.add_bonds(bonds)
    u.atoms.positions = np.array(xyz)
    return u


canvas = make_canvas()
helix = mn.Molecule(ideal_helix())
context = context_material(0.0)
helix.add_style("ball_and_stick", selection="not resid 5-9", material=context)
helix.add_style("ball_and_stick", selection="resid 5-9")
canvas.look_at(helix, viewpoint=(85, 0, 25), margin=0.15)

heading = title(helix, "Hydrogen bonds in an α helix", size=57)
subtitle = title(
    helix,
    "Carbonyl O(i) ··· H–N(i + 4)",
    position=(0.065, 0.145),
    size=30,
    color="#afbed0",
)
footer = title(
    helix,
    "Idealized backbone · 16 residues · explicit amide H · φ = −57° / ψ = −47°",
    position=(0.065, 0.93),
    size=21,
    color="#74869e",
)
count = title(
    helix, "12 expected|12 drawn", position=(0.065, 0.35), size=42, color=GOLD
)
# PM: HydrogenBonds(p, hydrogens="explicit", max_distance=3.5, min_angle=150).highlight(mode="3d")
network = []
for i in range(1, 13):
    ruler = helix.annotations.add_com_distance(
        selection1=f"resid {i} and name O",
        selection2=f"resid {i + 4} and name N",
        text1="",
        text2="",
    )
    ruler.visible = False
    network.append(ruler)
bond = helix.annotations.add_com_distance(
    selection1="resid 9 and name H",
    selection2="resid 5 and name O",
    text1="H",
    text2="O",
)
bond.visible = False
o_label = title(
    helix,
    "O · residue 5|Carbonyl acceptor · i",
    position=(0.72, 0.65),
    size=32,
    color="#f06478",
)
nh_label = title(
    helix,
    "N–H · residue 9|Amide donor · i + 4",
    position=(0.72, 0.32),
    size=32,
    color="#638cff",
)
turn = helix.get_view("resid 5-9")

with canvas.timeline() as t:
    show(t, heading, subtitle, footer)
    t.wait(1.8)
    show(t, *network, count)
    t.wait(2)
    t.play(canvas.camera.orbit(52), run_time=4)
    t.wait(0.6)
    show(t, *network, count, visible=False)
    t.play(
        t.tween(context.node.i.transparency, 0.94),  # SetOpacity(other, 0.06)
        canvas.camera.look_at(turn, margin=0.3),
        run_time=1.6,
    )
    show(t, bond, o_label, nh_label)
    t.play(canvas.camera.orbit(15), run_time=3)
    t.wait(1.2)
    show(t, bond, o_label, nh_label, visible=False)
    t.play(
        t.tween(context.node.i.transparency, 0.0),
        canvas.camera.look_at(helix, margin=0.15),
        run_time=1.6,
    )
    show(t, *network, count)
    t.play(canvas.camera.orbit(40), run_time=3)
    t.wait(0.8)

render(canvas, t, "alpha_helix_hbonds")
