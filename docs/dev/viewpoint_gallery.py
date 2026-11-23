"""
Render every automatic viewpoint method on ligand-bound structures from tests/data.

    uv run docs/dev/viewpoint_gallery.py OUT_DIR

Writes one PNG per structure and method, and results.json with the chosen
direction, solve time and how many ligand atoms have a clear ray to the camera
through the rendered cartoon. Pocket rows choose the direction from the ligand
and frame the ligand plus residues within 6 A; whole rows choose and frame on
the whole protein.
"""

import json
import os
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())
import bpy  # noqa: E402
import numpy as np  # noqa: E402
import molecularnodes as mn  # noqa: E402
from molecularnodes import viewpoint as vp  # noqa: E402
from molecularnodes.scene.base import _WorldBVH  # noqa: E402

DATA = Path(__file__).parents[2] / "tests" / "data"
CASES = {
    "1cbs": "REA",
    "1STP": "BTN",
    "1ATN": "ATP",
    "8H1B": "SAM",
    "4PNO": "U5P",
    "6ZPA": "3D1",
}
POCKET = ["front", "pca", "unobstructed", "visible_area", "raycast", "entropy"]
WHOLE = ["front", "pca", "entropy"]


def clear_rays(mol, points, forward):
    tree = _WorldBVH(mol.object, bpy.context.evaluated_depsgraph_get())
    return int(vp.ray_scorer([tree], points)(forward[None])[0])


def main(out: Path):
    out.mkdir(parents=True, exist_ok=True)
    canvas = mn.Canvas(resolution=(480, 360))
    canvas.engine = mn.scene.Cycles(samples=64, device="GPU")
    results = []
    for code, resname in CASES.items():
        for mode, methods in (("pocket", POCKET), ("whole", WHOLE)):
            for method in methods:
                mol = mn.Molecule.load(DATA / f"{code}.bcif")
                u = mol.universe
                ligand = u.select_atoms(f"resname {resname}").residues[0].atoms
                pocket = ligand | u.select_atoms(
                    "protein and around 6 group lig", lig=ligand
                )
                mol.add_style(
                    "cartoon", selection="protein", color=(0.6, 0.6, 0.65, 1.0)
                )
                mol.add_style("ball_and_stick", selection=ligand)
                lig_points = np.array(mol.get_view(ligand))
                target = lig_points if mode == "pocket" and method != "entropy" else mol
                start = time.perf_counter()
                canvas.look_at(target, viewpoint=method)
                seconds = time.perf_counter() - start
                canvas.look_at(mol.get_view(pocket) if mode == "pocket" else mol)
                forward = canvas.camera.basis[2]
                name = f"{code}_{mode}_{method}.png"
                canvas.snapshot(out / name)
                results.append(
                    {
                        "structure": code,
                        "ligand": resname,
                        "mode": mode,
                        "method": method,
                        "seconds": round(seconds, 3),
                        "clear": clear_rays(mol, lig_points, forward),
                        "atoms": len(lig_points),
                        "forward": [round(float(x), 3) for x in forward],
                        "image": name,
                    }
                )
                print(results[-1], flush=True)
                canvas.clear()
    (out / "results.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
