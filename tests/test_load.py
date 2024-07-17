import bpy
import numpy as np
import pytest
import itertools
import molecularnodes as mn
from .constants import data_dir, codes, attributes
from .utils import sample_attribute, NumpySnapshotExtension

mn._test_register()

styles = ["preset_1", "cartoon", "ribbon", "spheres", "surface", "ball_and_stick"]

centre_methods = ["", "centroid", "mass"]


def useful_function(snapshot_custom, style, code, assembly, cache_dir=None):
    obj = mn.entities.fetch(
        code, style=style, build_assembly=assembly, cache_dir=cache_dir
    ).object
    node = mn.blender.nodes.get_style_node(obj)

    if "Sphere As Points" in node.inputs.keys():
        node.inputs["Sphere As Points"].default_value = True

    mn.blender.nodes.realize_instances(obj)
    dont_realise = style == "cartoon" and code == "1BNA"
    for att in attributes:
        assert snapshot_custom == sample_attribute(obj, att, evaluate=dont_realise)


@pytest.mark.parametrize(
    "assembly, code, style", itertools.product([False], codes, styles)
)
def test_style_1(snapshot_custom: NumpySnapshotExtension, assembly, code, style):
    useful_function(snapshot_custom, style, code, assembly, cache_dir=data_dir)


# have to test a subset of styles with the biological assembly.
# testing some of the heavier styles run out of memory and fail on github actions


@pytest.mark.parametrize(
    "assembly, code, style",
    itertools.product([True], codes, ["cartoon", "surface", "ribbon"]),
)
def test_style_2(snapshot_custom: NumpySnapshotExtension, assembly, code, style):
    useful_function(snapshot_custom, style, code, assembly, cache_dir=data_dir)


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
        return mn.blender.mesh.named_attribute(object, "position")

    assert np.isclose(verts(mol), verts(mol2)).all()


@pytest.mark.parametrize("code, style", itertools.product(codes, styles))
def test_style_positions(snapshot_custom: NumpySnapshotExtension, code, style):
    mol = mn.entities.fetch(code, style=style, cache_dir=data_dir).object
    assert snapshot_custom == sample_attribute(mol, "position")


@pytest.mark.parametrize(
    "code, centre_method", itertools.product(codes, centre_methods)
)
def test_centring(snapshot_custom: NumpySnapshotExtension, code, centre_method):
    """fetch a pdb structure using code and translate the model using the
    centre_method. Check the CoG and CoM values against the snapshot file.
    """
    mol = mn.entities.fetch(code, centre=centre_method, cache_dir=data_dir)
    CoG = mol.centre()
    CoM = mol.centre(centre_type="mass")

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
        for method in centre_methods
    ]
    for mol1, mol2 in itertools.combinations(mols, 2):
        assert not np.allclose(
            mol1.centre(centre_type="centroid"), mol2.centre(centre_type="centroid")
        )
        assert not np.allclose(
            mol1.centre(centre_type="mass"), mol2.centre(centre_type="mass")
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
    for att in ["position"]:
        for mol in molecules:
            assert snapshot_custom == sample_attribute(mol, att, evaluate=False)


@pytest.mark.parametrize("del_hydrogen", [True, False])
def test_rcsb_nmr(snapshot_custom, del_hydrogen):
    mol = mn.entities.fetch(
        "2M6Q", style="cartoon", cache_dir=data_dir, del_hydrogen=del_hydrogen
    )
    assert len(mol.frames.objects) == 10
    assert (
        mol.object.modifiers["MolecularNodes"]
        .node_group.nodes["Animate Value"]
        .inputs["Value Max"]
        .default_value
        == 9
    )
    assert snapshot_custom == mol.named_attribute("position")
    assert snapshot_custom == sample_attribute(mol, "position", evaluate=True)

    bpy.context.scene.frame_set(1)
    pos_1 = mol.named_attribute("position", evaluate=True)
    bpy.context.scene.frame_set(100)
    pos_2 = mol.named_attribute("position", evaluate=True)
    bpy.context.scene.frame_set(1)
    assert (pos_1 != pos_2).all()


def test_load_small_mol(snapshot_custom):
    mol = mn.entities.load_local(data_dir / "ASN.cif")
    for att in ["position", "bond_type"]:
        assert snapshot_custom == sample_attribute(mol, att).tolist()


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
        obj_1 = mn.entities.fetch("6BQN", style="cartoon", cache_dir=test_cache)
        file = os.path.join(test_cache, "6BQN.bcif")
        assert os.path.exists(file)

        obj_2 = mn.entities.fetch("6BQN", style="cartoon", cache_dir=test_cache)
        assert (
            sample_attribute(obj_1, "position") == sample_attribute(obj_2, "position")
        ).all()
