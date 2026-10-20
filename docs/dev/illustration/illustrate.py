"""Render small-molecule vector illustrations with the prototype Style Illustration group.

uv run docs/dev/illustration/illustrate.py <outdir> [LIGAND ...]
"""

import os
import sys
import tempfile

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())
import bpy  # noqa: E402
import numpy as np  # noqa: E402
from biotite.structure.info import residue  # noqa: E402
from nodebpy import TreeBuilder  # noqa: E402
from nodebpy import geometry as g  # noqa: E402
import molecularnodes as mn  # noqa: E402
from molecularnodes import color as mncolor  # noqa: E402
from molecularnodes.converters.biotite import universe_from_atoms  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style_illustration import StyleIllustration  # noqa: E402

KEEP = "cpk"
LIGANDS = ["CFF", "HEM", "ATP", "NAG", "STI", "RET"]

# kwargs for StyleIllustration, plus "carbon" to draw carbon circles or not (skeletal look)
STYLES = {
    "classic": dict(fill=True, fill_element=0.0, atom_scale=0.28),
    "cpk": dict(
        fill=True,
        fill_element=1.0,
        stroke_color=(0.1, 0.1, 0.1, 1),
        bond_element=1.0,
        bond_thickness=0.18,
        atom_scale=0.35,
    ),
    "skeletal": dict(
        fill=True,
        carbon=False,
        fill_element=0.0,
        stroke_element=1.0,
        bond_element=1.0,
        atom_scale=0.22,
        bond_gap=0.05,
        outline_thickness=0.08,
    ),
    "outline": dict(
        fill=False,
        stroke_element=1.0,
        bond_element=1.0,
        atom_scale=0.3,
        bond_gap=0.1,
        bond_orders=False,
    ),
    "pastel": dict(
        fill=True,
        fill_element=0.45,
        stroke_element=1.0,
        bond_color=(0.35, 0.35, 0.4, 1),
        atom_scale=0.42,
        outline_thickness=0.1,
        bond_thickness=0.2,
        bond_orders=False,
    ),
}


def gp_material(name):
    mat = bpy.data.materials.get(name) or bpy.data.materials.new(name)
    if mat.grease_pencil is None:
        bpy.data.materials.create_gpencil_data(mat)
    mat.grease_pencil.color = (0, 0, 0, 1)
    mat.grease_pencil.fill_color = (1, 1, 1, 1)
    return mat


def ligand_array(code, hydrogens=False):
    arr = residue(code)
    if not hydrogens:
        arr = arr[arr.element != "H"]
    # align principal axes to x, y, z so planar groups lie in the xy plane
    xyz = arr.coord - arr.coord.mean(axis=0)
    _, _, vt = np.linalg.svd(xyz, full_matrices=False)
    if np.linalg.det(vt) < 0:
        vt[2] *= -1
    arr.coord = xyz @ vt.T
    arr.chain_id[:] = "A"
    return arr


def load(code):
    arr = ligand_array(code)
    mol = mn.Molecule(universe_from_atoms(arr), name=code)
    nums = mol.named_attribute("atomic_number")
    rgba = mncolor.colors_from_elements(nums) / 255
    rgba[nums == 6, :3] = 0.35
    mol["Color"] = rgba.astype(np.float32)
    with mol.tree.reset() as (atoms, join):
        atoms >> join
    return mol


def illustration(mol, camera, style):
    style = dict(style)
    carbon = style.pop("carbon", True)
    name = f"{mol.name}_illustration"
    with TreeBuilder.geometry(name, is_modifier=True, clear=True) as tree:
        out = tree.outputs.geometry("Grease Pencil")
        atoms = g.ObjectInfo(mol.object, transform_space="RELATIVE").o.geometry
        circles = (
            True
            if carbon
            else g.Compare.integer.not_equal(
                g.NamedAttribute.integer("atomic_number").o.attribute, 6
            )
        )
        (
            StyleIllustration(
                atoms,
                camera=camera,
                circles=circles,
                atom_material=gp_material("MN Atom"),
                bond_material=gp_material("MN Bond"),
                **style,
            )
            >> out
        )
    ob = bpy.data.objects.get(name)
    if ob is None:
        ob = bpy.data.objects.new(name, bpy.data.grease_pencils.new(name))
        bpy.context.scene.collection.objects.link(ob)
        ob.modifiers.new("Illustration", "NODES")
    ob.modifiers["Illustration"].node_group = bpy.data.node_groups[name]
    ob.data.materials.clear()
    for m in bpy.data.materials:
        if m.grease_pencil is not None:
            ob.data.materials.append(m)
    return ob


def main(out, ligands):
    canvas = mn.Canvas(mn.scene.EEVEE(samples=16), resolution=(500, 500))
    canvas.compositor.device = "CPU"
    canvas.view_transform = "Standard"
    canvas.background = (1, 1, 1, 1)
    camera = bpy.context.scene.camera
    for code in ligands:
        canvas.clear()
        mol = load(code)
        mol.object.hide_render = True
        canvas.camera.lens = 120
        canvas.look_at(mol, viewpoint="top", margin=0.15)
        for sname, style in STYLES.items():
            ob = illustration(mol, camera, style)
            canvas.snapshot(os.path.join(out, f"{code}_{sname}.png"))
            if sname == KEEP:
                bpy.ops.wm.save_as_mainfile(
                    filepath=os.path.join(out, f"{code}_{sname}.blend"), copy=True
                )
            bpy.data.objects.remove(ob)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:] or LIGANDS)
