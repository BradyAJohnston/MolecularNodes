import itertools
import bpy
import databpy as db
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.nodes.geometry import (
    StyleBallAndStick,
    StyleCartoon,
    StyleRibbon,
    StyleSpheres,
    StyleSurface,
)
from .constants import codes, data_dir
from .utils import GeometrySet, NumpySnapshotExtension

STYLES_TO_TEST = [
    "cartoon",
    "ribbon",
    "spheres",
    "surface",
    "ball_and_stick",
]

# `surface` builds its mesh through a volume, and the marching cubes step lands on a
# slightly different vertex count on each platform's Blender build. It is also by far
# the slowest style, so it is only exercised on the smallest structure.
SURFACE_CODE = "1BNA"
STYLE_PARAMS = [
    (code, assembly, style)
    for code, assembly, style in itertools.product(codes, [True, False], STYLES_TO_TEST)
    if style != "surface" or code == SURFACE_CODE
]


@pytest.mark.parametrize("code, assembly, style", STYLE_PARAMS)
def test_style_1(snapshot, code, assembly, style):
    mol = mn.Molecule.fetch(code, cache=data_dir)
    with mol.tree.reset() as (atoms, join):
        match style:
            case "ball_and_stick":
                style_node = StyleBallAndStick(sphere_geometry="Mesh")
            case "spheres":
                style_node = StyleSpheres(sphere_geometry="Mesh")
            case "cartoon":
                style_node = StyleCartoon()
            case "ribbon":
                style_node = StyleRibbon()
            case "surface":
                style_node = StyleSurface()

        assembly = (
            mn.nodes.geometry.AssemblyInstance(data_object=mol.create_data_object())
            if assembly
            else None
        )
        (atoms >> style_node >> assembly >> join)

    assert snapshot == GeometrySet(mol.object).summary()


@pytest.mark.parametrize(
    "code, format", list(itertools.product(codes, ["bcif", "cif", "pdb"]))
)
def test_download_format(code, format):
    mol = mn.Molecule.fetch(code, format=format, cache=data_dir)
    assert mol.props.entity_type == mn.entities.base.EntityType.MOLECULE.value
    with db.ObjectTracker() as o:
        bpy.ops.mn.import_molecule(
            code=code, file_format=format, cache_dir=str(data_dir)
        )
        mol2 = bpy.context.scene.MNSession.match(o.latest())

    assert np.allclose(mol.position, mol2.position)


@pytest.mark.parametrize("code", codes)
def test_style_positions(snapshot_custom: NumpySnapshotExtension, code):
    mol = mn.Molecule.fetch(code, cache=data_dir)
    assert snapshot_custom == mol.position


def test_local_pdb(snapshot_custom):
    molecules = [mn.Molecule.load(data_dir / f"1l58.{ext}") for ext in ("cif", "pdb")]
    molecules.append(mn.Molecule.fetch("1l58", format="bcif"))
    for mol in molecules:
        assert snapshot_custom == mol.named_attribute("position")


def test_pdb_blank_res_id(snapshot):
    # VESTA writes PDB files with blank res_id columns, which biotite refuses
    # to parse without help (#1091)
    mol = mn.Molecule.load(data_dir / "vesta_blank_res_id.pdb")
    assert mol.universe.atoms.n_atoms == 20
    assert (mol.universe.atoms.resids == 1).all()
    assert snapshot == mol.position


@pytest.mark.filterwarnings("ignore:.*elements were guessed.*:UserWarning")
def test_pdb_no_bonds(snapshot):
    mol = mn.Molecule.load(data_dir / "no_bonds.pdb")
    assert len(mol.object.data.edges) == 0
    assert snapshot == mol.position


def test_rcsb_nmr(snapshot_custom):
    mol = mn.Molecule.fetch("2M6Q", cache=data_dir)
    # multi-model (NMR) structures now load as trajectory frames of the Universe
    assert mol.universe.trajectory.n_frames > 1

    mol.add_style("cartoon")
    assert snapshot_custom == mol.named_attribute("position")

    # advancing the frame pushes the corresponding model's positions onto the mesh
    mol.set_frame(0)
    pos_0 = mol.named_attribute("position")
    mol.set_frame(5)
    pos_5 = mol.named_attribute("position")
    assert not np.allclose(pos_0, pos_5)


def test_load_small_mol(snapshot):
    mol = mn.Molecule.load(data_dir / "ASN.cif")
    assert mol.props.entity_type == mn.entities.base.EntityType.MOLECULE.value
    assert snapshot == GeometrySet(mol.object, strict=True).summary()
