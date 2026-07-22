import itertools
import bpy
import databpy as db
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.nodes.geometry import (
    AnimateFrames,
    AnimateValue,
    StyleBallAndStick,
    StyleCartoon,
    StyleRibbon,
    StyleSpheres,
    StyleSurface,
)
from .constants import attributes, codes, data_dir
from .utils import NumpySnapshotExtension

STYLES_TO_TEST = [
    "cartoon",
    "ribbon",
    "spheres",
    "surface",
    "ball_and_stick",
]


@pytest.mark.parametrize(
    "code, assembly, style", itertools.product(codes, [True, False], STYLES_TO_TEST)
)
def test_style_1(snapshot_custom: NumpySnapshotExtension, code, assembly, style):
    mol = mn.Molecule.fetch(code, cache=data_dir)
    with mol.tree.reset() as (atoms, join):
        match style:
            case "ball_and_stick":
                style_node = StyleBallAndStick(sphere_geometry="Mesh")
            case "spheres":
                style_node = StyleSpheres(geometry="Mesh")
            case "cartoon":
                style_node = StyleCartoon()
            case "ribbon":
                style_node = StyleRibbon()
            case "surface":
                style_node = StyleSurface()

        assembly = (
            mn.nodes.geometry.AssemblyInstance(
                data_object=mol.create_data_object(), realize_all=True
            )
            if assembly
            else None
        )
        (
            atoms
            >> style_node
            >> assembly
            >> join
        )

    for att in attributes:
        try:
            assert snapshot_custom == mol.named_attribute(
                att, evaluate=style == "cartoon" and code == "1BNA"
            )
        except AttributeError as e:
            assert snapshot_custom == e


@pytest.mark.parametrize(
    "code, format", itertools.product(codes, ["bcif", "cif", "pdb"])
)
def test_download_format(code, format):
    mol = mn.Molecule.fetch(code, format=format, cache=data_dir)
    assert mol._entity_type == mn.entities.base.EntityType.MOLECULE
    assert mol.object.mn.entity_type == mol._entity_type.value
    with db.ObjectTracker() as o:
        bpy.ops.mn.import_fetch(code=code, file_format=format, cache_dir=str(data_dir))
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


@pytest.mark.filterwarnings("ignore:.*elements were guessed.*:UserWarning")
def test_pdb_no_bonds(snapshot):
    mol = mn.Molecule.load(data_dir / "no_bonds.pdb")
    assert len(mol.object.data.edges) == 0
    assert snapshot == mol.position


def test_rcsb_nmr(snapshot_custom):
    mol = mn.Molecule.fetch("2M6Q", cache=data_dir)

    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> AnimateFrames(frames=mol.frames, frame=AnimateValue(value_max=9))
            >> StyleCartoon()
            >> join
        )
    assert snapshot_custom == mol.named_attribute("position")

    bpy.context.scene.frame_set(1)
    pos_1 = mol.named_attribute("position", evaluate=True)
    bpy.context.scene.frame_set(100)
    pos_2 = mol.named_attribute("position", evaluate=True)
    bpy.context.scene.frame_set(1)
    assert not np.allclose(pos_1, pos_2)


def test_load_small_mol(snapshot_custom):
    mol = mn.Molecule.load(data_dir / "ASN.cif")
    assert mol._entity_type == mn.entities.base.EntityType.MOLECULE
    assert mol.object.mn.entity_type == mol._entity_type.value
    for att in ["position", "bond_type"]:
        assert snapshot_custom == mol.named_attribute(att).tolist()
