import itertools
import bpy
import databpy as db
import numpy as np
import pytest
import molecularnodes as mn
from .constants import attributes, codes, data_dir
from .utils import NumpySnapshotExtension

STYLES_TO_TEST = [
    "preset_1",
    "cartoon",
    "ribbon",
    "spheres",
    "surface",
    "ball_and_stick",
]

CENTRE_METHODS_TO_TEST = ["", "centroid", "mass"]


@pytest.mark.parametrize(
    "assembly, code, style", itertools.product([False], codes, STYLES_TO_TEST)
)
def test_style_1(snapshot_custom: NumpySnapshotExtension, assembly, code, style):
    mol = mn.Molecule.fetch(code, cache=data_dir).add_style(
        style=style, assembly=assembly
    )
    if style == "spheres":
        mol.styles[0].geometry = "Mesh"

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


@pytest.mark.parametrize(
    "code, centre_method", itertools.product(codes, CENTRE_METHODS_TO_TEST)
)
def test_centring(snapshot_custom: NumpySnapshotExtension, code, centre_method):
    """fetch a pdb structure using code and translate the model using the
    centre_method. Check the CoG and CoM values against the snapshot file.
    """
    mol = mn.Molecule.fetch(code, cache=data_dir).centre_molecule(centre_method)
    CoG = mol.centroid()
    CoM = mol.centroid(weight="mass")

    if centre_method == "centroid":
        assert np.linalg.norm(CoG) < 1e-06
    elif centre_method == "mass":
        assert np.linalg.norm(CoM) < 1e-06

    CoG = np.array_str(CoG, precision=4, suppress_small=True)
    CoM = np.array_str(CoM, precision=4, suppress_small=True)
    assert snapshot_custom == [CoG, CoM]


@pytest.mark.parametrize("code", codes)
def test_centring_different(code):
    """fetch multiple instances of the same pdb structure and translate
    each by a different centring method. Check that their centroids and
    positions are in fact different.
    """
    mols = [
        mn.Molecule.fetch(code, cache=data_dir).centre_molecule(method)
        for method in CENTRE_METHODS_TO_TEST
    ]
    for mol1, mol2 in itertools.combinations(mols, 2):
        assert not np.allclose(mol1.centroid(), mol2.centroid())
        assert not np.allclose(
            mol1.centroid(weight="mass"), mol2.centroid(weight="mass")
        )
        assert not np.allclose(
            mol1.named_attribute("position"), mol2.named_attribute("position")
        )


def test_local_pdb(snapshot_custom):
    molecules = [mn.Molecule.load(data_dir / f"1l58.{ext}") for ext in ("cif", "pdb")]
    molecules.append(mn.Molecule.fetch("1l58", format="bcif"))
    for mol in molecules:
        assert snapshot_custom == mol.named_attribute("position")


def test_pdb_no_bonds(snapshot):
    mol = mn.Molecule.load(data_dir / "no_bonds.pdb")
    assert len(mol.object.data.edges) == 0
    assert snapshot == mol.position


def test_rcsb_nmr(snapshot_custom):
    mol = mn.Molecule.fetch("2M6Q", cache=data_dir).add_style(style="cartoon")
    assert len(mol.frames.objects) == 10
    assert mol.node_group.nodes["Animate Value"].inputs["Value Max"].default_value == 9
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
