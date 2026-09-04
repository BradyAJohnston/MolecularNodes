import itertools
import shutil
import bpy
import numpy as np
import pytest
from databpy import ObjectTracker
import molecularnodes as mn
from .constants import codes, data_dir
from .utils import NumpySnapshotExtension


@pytest.mark.parametrize("code", codes)
def test_op_fetch(snapshot_custom: NumpySnapshotExtension, code):
    scene = bpy.context.scene
    style = "ribbon"
    format = "cif"

    with ObjectTracker() as o:
        bpy.ops.mn.import_molecule(code=code, file_format=format, style=style)
        mol1 = scene.MNSession.match(o.latest())

    with ObjectTracker() as o:
        bpy.ops.mn.import_molecule(
            code=f"pdb_{code.rjust(8, '0')}", file_format=format, style=style
        )
        mol2 = scene.MNSession.match(o.latest())

    mol3 = mn.Molecule.fetch(code, format=format, cache=data_dir)
    mol3.add_style(style=style)

    for test1, test2 in itertools.combinations([mol1, mol2, mol3], 2):
        np.testing.assert_allclose(test1.position, test2.position)


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("file_format", ["bcif", "cif", "pdb"])
def test_op_local(snapshot_custom, code, file_format):
    session = bpy.context.scene.MNSession
    path = mn.download.StructureDownloader(cache=data_dir).download(
        code=code, format=file_format
    )

    with ObjectTracker() as o:
        bpy.ops.mn.import_molecule(method="local", filepath=str(path), node_setup=False)
        mol = session.match(o.latest())
        assert mol.props.entity_type == mn.entities.base.EntityType.MOLECULE.value

    assert snapshot_custom == mol.position


def test_op_dropped_files():
    """The file handler invokes mn.import_molecule with directory + files set."""
    session = bpy.context.scene.MNSession
    paths = [
        mn.download.StructureDownloader(cache=data_dir).download(
            code=code, format="cif"
        )
        for code in codes[:2]
    ]

    with ObjectTracker() as o:
        res = bpy.ops.mn.import_molecule(
            directory=str(paths[0].parent),
            files=[{"name": path.name} for path in paths],
            style="ribbon",
        )
        assert res == {"FINISHED"}
        objects = o.new_objects()

    mols = [session.match(obj) for obj in objects]
    assert len(mols) == len(paths)
    for mol, path in zip(mols, paths):
        assert mol is not None
        assert mol.name == path.name
        # dropped files get the same node setup as dialog imports
        style_trees = [
            node.node_tree.name
            for node in mol.modifier_node_tree.nodes
            if isinstance(node, bpy.types.GeometryNodeGroup)
        ]
        assert "Style Ribbon" in style_trees
        assert "Set Color" in style_trees


def _build_structure_folder(tmp_path, name):
    """A folder of structure files with nested subfolders, plus files to skip."""
    src = data_dir / "1BNA.pdb"
    root = tmp_path / name
    (root / "groupA" / "deep").mkdir(parents=True)
    (root / "groupB").mkdir()
    (root / "empty").mkdir()
    shutil.copy(src, root / "top.pdb")
    shutil.copy(src, root / "groupA" / "a1.pdb")
    shutil.copy(src, root / "groupA" / "deep" / "d1.pdb")
    shutil.copy(src, root / "groupB" / "b1.pdb")
    (root / "groupA" / "notes.txt").write_text("not a structure")
    return root


def test_op_import_folder(tmp_path):
    """A folder filepath imports every structure file in it recursively, placing
    the objects in nested collections that mirror the folder layout."""
    session = bpy.context.scene.MNSession
    root = _build_structure_folder(tmp_path, "folder_import")

    with ObjectTracker() as o:
        res = bpy.ops.mn.import_molecule(
            method="local", filepath=str(root), node_setup=False
        )
        assert res == {"FINISHED"}
        objects = o.new_objects()

    assert len(objects) == 4
    for obj in objects:
        assert session.match(obj) is not None

    root_coll = bpy.data.collections["folder_import"]
    assert root_coll.name in bpy.data.collections["Molecular Nodes"].children
    group_a = root_coll.children["groupA"]
    assert "top.pdb" in root_coll.objects
    assert "a1.pdb" in group_a.objects
    assert "d1.pdb" in group_a.children["deep"].objects
    assert "b1.pdb" in root_coll.children["groupB"].objects
    # folders without structure files get no collection, non-structure files skip
    assert "empty" not in root_coll.children

    # a second import into the same folder reuses the collections
    with ObjectTracker() as o:
        bpy.ops.mn.import_molecule(method="local", filepath=str(root), node_setup=False)
        assert len(o.new_objects()) == 4
    assert "folder_import.001" not in bpy.data.collections
    assert len(group_a.objects) == 2


def test_op_import_folder_objects_only(tmp_path):
    """With objects_only set, the folder import keeps the created objects but
    discards the Python entities, which stay individually relinkable."""
    session = bpy.context.scene.MNSession
    root = _build_structure_folder(tmp_path, "folder_objects_only")

    with ObjectTracker() as o:
        res = bpy.ops.mn.import_molecule(
            method="local", filepath=str(root), node_setup=False, objects_only=True
        )
        assert res == {"FINISHED"}
        objects = o.new_objects()

    assert len(objects) == 4
    for obj in objects:
        assert session.match(obj) is None
        assert mn.entities.reload.can_reload(obj)

    # a discarded entity can be reconstructed from the object's recorded source
    entity = mn.entities.reload.reload_entity(objects[0])
    assert session.match(objects[0]) is entity


def test_op_api_mda(snapshot_custom: NumpySnapshotExtension):
    bpy.context.scene.frame_set(0)

    topo = str(data_dir / "md_ppr/box.gro")
    traj = str(data_dir / "md_ppr/first_5_frames.xtc")
    name = "AnotherNewTrajectory"

    with ObjectTracker() as o:
        bpy.ops.mn.import_molecule(
            method="local",
            filepath=topo,
            trajectory=traj,
            style="ribbon",
        )
        obj_1 = o.latest()

    traj_op = bpy.context.scene.MNSession.match(obj_1)
    assert traj_op.name == "box.gro"
    traj_op.name = name
    assert traj_op.name == name
    assert traj_op._mn_entity_type == mn.entities.base.EntityType.MOLECULE

    traj_func = mn.entities.molecule.Molecule.load(
        topo, traj, name="test", style="ribbon"
    )

    bpy.context.scene.frame_set(2)
    assert np.allclose(traj_func.position, traj_op.position)
    pos_2 = traj_func.position.copy()
    bpy.context.scene.frame_set(4)
    traj_op.set_frame(4)
    traj_func.set_frame(4)

    assert not np.allclose(pos_2, traj_op.position)
    assert not np.allclose(pos_2, traj_func.position)


def test_op_dropped_files_share_node_group():
    """With Share Node Group enabled, all dropped structures use one tree."""
    session = bpy.context.scene.MNSession
    paths = [
        mn.download.StructureDownloader(cache=data_dir).download(
            code=code, format="cif"
        )
        for code in codes[:2]
    ]

    with ObjectTracker() as o:
        res = bpy.ops.mn.import_molecule(
            directory=str(paths[0].parent),
            files=[{"name": path.name} for path in paths],
            style="ribbon",
            share_node_group=True,
        )
        assert res == {"FINISHED"}
        objects = o.new_objects()

    mols = [session.match(obj) for obj in objects]
    trees = [mol.modifier_node_tree for mol in mols]
    assert trees[0] == trees[1]
    style_trees = [
        node.node_tree.name
        for node in trees[0].nodes
        if isinstance(node, bpy.types.GeometryNodeGroup)
    ]
    assert "Style Ribbon" in style_trees


def test_share_node_group_api():
    """Assigning one molecule's tree to another shares styling between them."""
    mol1 = mn.Molecule.fetch(codes[0], cache=data_dir, format="cif")
    mol2 = mn.Molecule.fetch(codes[1], cache=data_dir, format="cif")
    mol1.add_style("spheres")

    mol2.tree = mol1.tree
    assert mol2.modifier_node_tree == mol1.modifier_node_tree

    # a style added through either entity lands in the shared tree
    mol2.add_style("ribbon")
    style_trees = [
        node.node_tree.name
        for node in mol1.modifier_node_tree.nodes
        if isinstance(node, bpy.types.GeometryNodeGroup)
    ]
    assert "Style Spheres" in style_trees
    assert "Style Ribbon" in style_trees


def test_op_fetch_comma_separated_codes():
    """Comma-separated codes fetch multiple structures, sharing one tree."""
    session = bpy.context.scene.MNSession

    with ObjectTracker() as o:
        res = bpy.ops.mn.import_molecule(
            code=f"{codes[0]}, {codes[1]}",
            file_format="cif",
            style="ribbon",
            share_node_group=True,
            cache_dir=str(data_dir),
        )
        assert res == {"FINISHED"}
        objects = o.new_objects()

    assert len(objects) == 2
    mols = [session.match(obj) for obj in objects]
    assert sorted(mol.props.code for mol in mols) == sorted(codes[:2])
    trees = [mol.modifier_node_tree for mol in mols]
    assert trees[0] == trees[1]

    # without sharing, each fetched structure gets its own tree
    with ObjectTracker() as o:
        res = bpy.ops.mn.import_molecule(
            code=f"{codes[0]},{codes[1]}",
            file_format="cif",
            style="ribbon",
            cache_dir=str(data_dir),
        )
        assert res == {"FINISHED"}
        objects = o.new_objects()

    mols = [session.match(obj) for obj in objects]
    trees = [mol.modifier_node_tree for mol in mols]
    assert trees[0] != trees[1]
