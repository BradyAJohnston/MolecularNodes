import bpy
import numpy as np
import pytest
from databpy import ObjectTracker
import molecularnodes as mn
from .constants import attributes, codes, data_dir
from .utils import NumpySnapshotExtension


@pytest.mark.parametrize("code", codes)
def test_op_api_cartoon(
    snapshot_custom: NumpySnapshotExtension, code, style="ribbon", format="bcif"
):
    scene = bpy.context.scene

    with ObjectTracker() as o:
        bpy.ops.mn.import_fetch(code=code, file_format=format, style=style)
        mol1 = scene.MNSession.match(o.latest())

    mol2 = mn.Molecule.fetch(code, format=format, cache=data_dir)
    mol2.add_style(style=style)

    # objects being imported via each method should have identical snapshots
    for mol in [mol1, mol2]:
        for name in attributes:
            try:
                assert snapshot_custom == mol.named_attribute(
                    "position", evaluate=False
                )
            except AttributeError as e:
                assert snapshot_custom == str(e)


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("file_format", ["bcif", "cif", "pdb"])
def test_op_local(snapshot_custom, code, file_format):
    session = bpy.context.scene.MNSession
    path = mn.download.download(code=code, format=file_format, cache=data_dir)

    with ObjectTracker() as o:
        bpy.ops.mn.import_local(filepath=str(path), node_setup=False)
        mol = session.match(o.latest())

    with ObjectTracker() as o:
        bpy.ops.mn.import_local(filepath=path, centre=True, centre_type="centroid")
        mol_cent = session.match(o.latest())

    assert snapshot_custom == mol.position
    assert snapshot_custom == mol_cent.position
    assert not np.allclose(mol.position, mol_cent.position)


def test_op_api_mda(snapshot_custom: NumpySnapshotExtension):
    bpy.context.scene.frame_set(0)

    topo = str(data_dir / "md_ppr/box.gro")
    traj = str(data_dir / "md_ppr/first_5_frames.xtc")
    name = "AnotherNewTrajectory"

    with ObjectTracker() as o:
        bpy.ops.mn.import_trajectory(
            topology=topo, trajectory=traj, name=name, style="ribbon"
        )
        obj_1 = o.latest()

    traj_op = bpy.context.scene.MNSession.match(obj_1)
    assert traj_op.name == name

    traj_func = mn.entities.trajectory.load(topo, traj, name="test", style="ribbon")

    bpy.context.scene.frame_set(2)
    assert np.allclose(traj_func.position, traj_op.position)
    pos_2 = traj_func.position.copy()
    bpy.context.scene.frame_set(4)
    traj_op.set_frame(4)
    traj_func.set_frame(4)

    assert not np.allclose(pos_2, traj_op.position)
    assert not np.allclose(pos_2, traj_func.position)


@pytest.mark.skipif(
    bpy.app.version_string.startswith("4.2"),
    reason="Test fails in 4.2 but succeeeds otherwise",
)
def test_op_residues_selection_custom():
    topo = str(data_dir / "md_ppr/box.gro")
    traj = str(data_dir / "md_ppr/first_5_frames.xtc")

    with ObjectTracker() as o:
        bpy.ops.mn.import_trajectory(
            topology=topo, trajectory=traj, name="NewTrajectory", style="ribbon"
        )
    area = bpy.context.screen.areas[-1]
    area.ui_type = "GeometryNodeTree"
    with bpy.context.temp_override(area=area):
        bpy.ops.mn.residues_selection_custom("EXEC_DEFAULT")
