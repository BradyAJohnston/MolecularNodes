import bpy
import os
import pytest
import molecularnodes as mn
from molecularnodes.blender.mesh import named_attribute

import MDAnalysis as mda
import numpy as np
from .constants import data_dir
from .utils import sample_attribute, NumpySnapshotExtension

mn._test_register()


class TestTrajectory:
    @pytest.fixture(scope="module")
    def universe(self):
        top = data_dir / "md_ppr/box.gro"
        traj = data_dir / "md_ppr/first_5_frames.xtc"
        u = mda.Universe(top, traj)
        return u

    @pytest.fixture(scope="module")
    def universe_with_bonds(self):
        top = data_dir / "md_ppr/md.tpr"
        traj = data_dir / "md_ppr/md.gro"
        u = mda.Universe(top, traj)
        return u

    @pytest.fixture(scope="module")
    def Trajectory_cross_boundary(self):
        topo = data_dir / "martini/dode_membrane/topol_nowat.gro"
        traj = data_dir / "martini/dode_membrane/traj_imaged_dt1ns_frames_1-10.xtc"
        u = mda.Universe(topo, traj)
        traj = mn.entities.Trajectory(u)
        traj.create_object()
        return traj

    @pytest.fixture(scope="module")
    def Trajectory(self, universe):
        traj = mn.entities.Trajectory(universe)
        traj.create_object()
        return traj

    @pytest.fixture(scope="module")
    def Trajectory_with_bonds(self, universe_with_bonds):
        traj = mn.entities.Trajectory(universe_with_bonds)
        traj.create_object()
        return traj

    @pytest.fixture(scope="module")
    def session(self):
        return mn.session.get_session()

    def test_include_bonds(self, Trajectory_with_bonds):
        assert Trajectory_with_bonds.object.data.edges.items() != []

    def test_attributes_added(self, Trajectory):
        attributes = Trajectory.object.data.attributes.keys()
        # check if all attributes are added.

        attribute_added = [
            "vdw_radii",
            "b_factor",
            "atomic_number",
            "res_id",
            "res_name",
            "chain_id",
            "atom_types",
            "atom_name",
            "position",
            "is_backbone",
            "is_alpha_carbon",
            "is_solvent",
            "is_nucleic",
            "is_peptide",
        ]
        for att in attribute_added:
            assert att in attributes

    def test_trajectory_update(self, snapshot_custom, Trajectory):
        traj = Trajectory
        bpy.context.scene.frame_set(0)
        pos_a = traj.named_attribute("position")
        assert snapshot_custom == pos_a

        bpy.context.scene.frame_set(4)
        pos_b = traj.named_attribute("position")
        assert snapshot_custom == pos_b

        assert not np.allclose(pos_a, pos_b)

    @pytest.mark.parametrize("offset", [-2, 2])
    def test_trajectory_offset(self, Trajectory, offset):
        traj = Trajectory
        traj.offset = 0
        bpy.context.scene.frame_set(0)
        pos_0 = traj.named_attribute("position")

        # if the offset is negative, the positions of the starting frame 0 will change.
        # if the offset is positive, then all of the frames up till the offset frame
        # will remain the same
        traj.offset = offset
        if offset < 0:
            assert not np.allclose(pos_0, traj.named_attribute("position"))
        else:
            assert np.allclose(pos_0, traj.named_attribute("position"))
            bpy.context.scene.frame_set(offset - 1)
            assert np.allclose(pos_0, traj.named_attribute("position"))

        # after resetting the offset to 0, it should be the same as the initial positions
        bpy.context.scene.frame_set(0)
        traj.offset = 0
        assert np.allclose(pos_0, traj.named_attribute("position"))

    @pytest.mark.parametrize("interpolate", [True, False])
    def test_subframes(self, Trajectory, interpolate):
        traj = Trajectory
        bpy.context.scene.frame_set(0)
        traj.subframes = 0
        traj.interpolate = interpolate
        verts_a = traj.named_attribute("position")

        bpy.context.scene.frame_set(1)
        verts_b = traj.named_attribute("position")

        # should be different because we have changed the frame
        assert not np.allclose(verts_a, verts_b)

        for subframes in [1, 2, 3, 4]:
            bpy.context.scene.frame_set(0)
            frame = 1
            fraction = frame % (subframes + 1) / (subframes + 1)
            traj.subframes = subframes

            bpy.context.scene.frame_set(frame)
            verts_c = traj.named_attribute("position")

            if interpolate:
                # now using subframes and having interpolate=True there should be a difference
                assert not np.allclose(verts_b, verts_c)
                assert np.allclose(verts_c, mn.utils.lerp(verts_a, verts_b, t=fraction))
            else:
                # without using interopolation, the subframes means it should default back
                # to the previous best selected frame
                assert np.allclose(verts_a, verts_c)

    def test_correct_periodic(self, snapshot_custom, Trajectory_cross_boundary):
        u = Trajectory_cross_boundary
        u.subframes = 5
        bpy.context.scene.frame_set(2)
        pos_a = u.named_attribute("position")
        u.object.mn.correct_periodic = False
        pos_b = u.named_attribute("position")

        assert not np.allclose(pos_a, pos_b)
        assert snapshot_custom == pos_a

    def test_update_selection(self, snapshot_custom, Trajectory):
        # to API add selections we currently have to operate on the UIList rather than the
        # universe itself, which isn't great
        u = Trajectory
        bpy.context.scene.frame_set(0)
        sel = u.add_selection(name="custom_sel_1", selection_str="around 3.5 protein")
        bpy.context.scene.frame_set(5)
        sel_1 = u.named_attribute("custom_sel_1")
        bpy.context.scene.frame_set(50)
        sel_2 = u.named_attribute("custom_sel_1")
        # when we are updating, the selection around the protein will change from frame
        # to frame
        assert not (sel_1 != sel_2).all()

        # if we stop the selection from updating, then even when we change the frame
        # the selection will remain the same
        sel.updating = False
        bpy.context.scene.frame_set(100)
        assert (sel_2 == u.named_attribute("custom_sel_1")).all()
        # if we change the selection to updating, then the selection will be updated
        # and will no longer match with what came earlier
        sel.updating = False
        assert not (sel_2 != u.named_attribute("custom_sel_1")).all()

    def test_save_persistance(
        self,
        snapshot_custom: NumpySnapshotExtension,
        tmp_path,
        universe,
        session: mn.session.MNSession,
    ):
        session.clear()
        traj = mn.entities.Trajectory(universe)
        traj.create_object()
        uuid = traj.uuid
        bpy.context.scene.frame_set(0)
        filepath = str(tmp_path / "test.blend")

        # test that we can save the file and it is created only after saving
        assert not os.path.exists(session.stashpath(filepath))
        bpy.ops.wm.save_as_mainfile(filepath=filepath)
        assert os.path.exists(filepath)
        assert os.path.exists(session.stashpath(filepath))
        bpy.ops.wm.open_mainfile(filepath=filepath)

        traj = mn.session.get_session().trajectories[uuid]
        verts_frame_0 = traj.named_attribute("position")
        bpy.context.scene.frame_set(4)
        verts_frame_4 = traj.named_attribute("position")

        assert snapshot_custom == verts_frame_4
        assert not np.allclose(verts_frame_0, verts_frame_4)


@pytest.mark.parametrize("toplogy", ["pent/prot_ion.tpr", "pent/TOPOL2.pdb"])
def test_martini(snapshot_custom: NumpySnapshotExtension, toplogy):
    universe = mda.Universe(
        data_dir / "martini" / toplogy, data_dir / "martini/pent/PENT2_100frames.xtc"
    )
    traj = mn.entities.Trajectory(universe)
    traj.create_object()
    obj = traj.object
    bpy.context.scene.frame_set(0)
    pos_a = traj.named_attribute("position")

    bpy.context.scene.frame_set(50)
    pos_b = traj.named_attribute("position")
    assert not np.allclose(pos_a, pos_b)

    for att in obj.data.attributes.keys():
        assert snapshot_custom == traj.named_attribute(att)

    for att in obj.data.attributes.keys():
        assert snapshot_custom == traj.named_attribute(att)
