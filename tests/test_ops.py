import itertools
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
    assert traj_op._mn_entity_type == mn.entities.base.EntityType.MD.value

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
