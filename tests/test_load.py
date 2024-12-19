import bpy
import numpy as np
import pytest
import itertools
import molecularnodes as mn
import databpy as db
from .constants import data_dir, codes, attributes
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
    # have to test a subset of styles with the biological assembly.
    # testing some of the heavier styles run out of memory and fail on github actions
    if assembly:
        styles = ["cartoon", "surface", "ribbon"]

    mol = mn.entities.fetch(
        code, style=style, build_assembly=assembly, cache_dir=data_dir
    )
    node = mn.blender.nodes.get_style_node(mol.object)

    if "Sphere As Mesh" in node.inputs.keys():
        node.inputs["Sphere As Mesh"].default_value = True

    mn.blender.nodes.realize_instances(mol.object)
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
    mol = mn.entities.fetch(code, format=format, style=None, cache_dir=data_dir).object
    scene = bpy.context.scene
    scene.MN_pdb_code = code
    scene.MN_import_node_setup = False
    scene.MN_import_format_download = format
    names = [o.name for o in bpy.data.objects]
    bpy.ops.mn.import_wwpdb()

    for o in bpy.data.objects:
        if o.name not in names:
            mol2 = o

    def verts(object):
        return db.named_attribute(object, "position")

    assert np.isclose(verts(mol), verts(mol2)).all()


@pytest.mark.parametrize("code", codes)
def test_style_positions(snapshot_custom: NumpySnapshotExtension, code):
    mol = mn.entities.fetch(code, style=None, cache_dir=data_dir)
    assert snapshot_custom == mol.position


@pytest.mark.parametrize(
    "code, centre_method", itertools.product(codes, CENTRE_METHODS_TO_TEST)
)
def test_centring(snapshot_custom: NumpySnapshotExtension, code, centre_method):
    """fetch a pdb structure using code and translate the model using the
    centre_method. Check the CoG and CoM values against the snapshot file.
    """
    mol = mn.entities.fetch(code, centre=centre_method, cache_dir=data_dir)
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
        mn.entities.fetch(code, centre=method, cache_dir=data_dir)
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


# THESE TEST FUNCTIONS ARE NOT RUN
def test_local_pdb(snapshot_custom):
    molecules = [
        mn.entities.load_local(data_dir / f"1l58.{ext}", style="spheres")
        for ext in ("cif", "pdb")
    ]
    molecules.append(mn.entities.fetch("1l58", format="bcif"))
    for mol in molecules:
        assert snapshot_custom == mol.named_attribute("position")


def test_pdb_no_bonds(snapshot):
    mol = mn.entities.load_local(data_dir / "no_bonds.pdb", style=None)
    assert len(mol.object.data.edges) == 0
    assert snapshot == mol.position


@pytest.mark.parametrize("del_hydrogen", [True, False])
def test_rcsb_nmr(snapshot_custom, del_hydrogen):
    mol = mn.entities.fetch(
        "2M6Q", style="cartoon", cache_dir=data_dir, del_hydrogen=del_hydrogen
    )
    assert len(mol.frames.objects) == 10
    assert mol.node_group.nodes["Animate Value"].inputs["Value Max"].default_value == 9
    assert snapshot_custom == mol.named_attribute("position")

    bpy.context.scene.frame_set(1)
    pos_1 = mol.named_attribute("position", evaluate=True)
    bpy.context.scene.frame_set(100)
    pos_2 = mol.named_attribute("position", evaluate=True)
    bpy.context.scene.frame_set(1)
    assert (pos_1 != pos_2).all()


def test_load_small_mol(snapshot_custom):
    mol = mn.entities.load_local(data_dir / "ASN.cif")
    for att in ["position", "bond_type"]:
        assert snapshot_custom == mol.named_attribute(att).tolist()


def test_rcsb_cache(snapshot_custom):
    from pathlib import Path
    import tempfile
    import os

    # we want to make sure cached files are freshly downloaded, but
    # we don't want to delete our entire real cache
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as data_dir:
        test_cache = Path(data_dir)

        # Run the test
        mol1 = mn.entities.fetch("6BQN", style="cartoon", cache_dir=test_cache)
        file = os.path.join(test_cache, "6BQN.bcif")
        assert os.path.exists(file)

        mol2 = mn.entities.fetch("6BQN", style="cartoon", cache_dir=test_cache)
        assert np.allclose(mol1.position, mol2.position)
