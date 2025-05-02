import itertools
import os
import bpy
import databpy
import MDAnalysis as mda
import numpy as np
import pytest
import molecularnodes as mn
from .constants import data_dir
from .utils import NumpySnapshotExtension


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
    def univ_across_boundary(self):
        topo = data_dir / "martini/dode_membrane/topol_nowat.gro"
        traj = data_dir / "martini/dode_membrane/traj_imaged_dt1ns_frames_1-10.xtc"
        u = mda.Universe(topo, traj)
        return u

    @pytest.fixture(scope="module")
    def session(self):
        return mn.session.get_session()

    def test_include_bonds(self, universe_with_bonds):
        traj = mn.entities.Trajectory(universe_with_bonds)
        traj.create_object()
        assert traj.edges.items() != []

    def test_attributes_added(self, universe):
        traj = mn.entities.Trajectory(universe)
        traj.create_object()
        attributes = traj.list_attributes()
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

    def test_trajectory_update(self, snapshot, universe):
        traj = mn.entities.Trajectory(universe)
        traj.create_object(name="TestTrajectoryUpdate")
        print(f"{bpy.context.scene.frame_current=}")
        print(f"{list(bpy.app.handlers.frame_change_pre)=}")
        bpy.context.scene.frame_set(0)
        pos_a = traj.position
        assert snapshot == pos_a
        bpy.context.scene.frame_set(3)
        pos_b = traj.position
        print(f"{bpy.context.scene.MNSession.entities.keys()=}")
        [
            print("\n\n{}: {}".format(v._object_name, v.uuid))
            for v in bpy.context.scene.MNSession.entities.values()
        ]
        [print("{}: {}".format(obj, obj.uuid)) for obj in bpy.data.objects]
        print(f"{bpy.context.scene.MNSession.entities.values()=}")
        assert not np.allclose(pos_a, pos_b)
        assert snapshot == pos_b

    @pytest.mark.parametrize("offset", [-2, 2])
    def test_trajectory_offset(self, universe, offset: bool):
        traj = mn.entities.Trajectory(universe)
        traj.create_object()
        bpy.context.scene.frame_set(0)
        pos_0 = traj.position

        # if the offset is negative, the positions of the starting frame 0 will change.
        # if the offset is positive, then all of the frames up till the offset frame
        # will remain the same
        traj.offset = offset
        if offset < 0:
            assert not np.allclose(pos_0, traj.position)
        else:
            assert np.allclose(pos_0, traj.position)
            bpy.context.scene.frame_set(4)
            assert not np.allclose(pos_0, traj.position)

        # after resetting the offset to 0, it should be the same as the initial positions
        bpy.context.scene.frame_set(0)
        traj.offset = 0
        assert np.allclose(pos_0, traj.position)

    @pytest.mark.parametrize("interpolate", [True, False])
    def test_subframes(self, universe, interpolate: bool):
        traj = mn.entities.Trajectory(universe)
        traj.create_object()
        bpy.context.scene.frame_set(0)
        traj.subframes = 0
        traj.interpolate = interpolate
        verts_a = traj.position

        bpy.context.scene.frame_set(1)
        verts_b = traj.position

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
                assert np.allclose(verts_c, databpy.lerp(verts_a, verts_b, t=fraction))
            else:
                # without using interopolation, the subframes means it should default back
                # to the previous best selected frame
                assert np.allclose(verts_a, verts_c)

    def test_correct_periodic(
        self,
        snapshot_custom: NumpySnapshotExtension,
        univ_across_boundary,
    ):
        traj = mn.entities.Trajectory(univ_across_boundary)
        traj.create_object()
        traj.subframes = 5
        bpy.context.scene.frame_set(2)
        pos_a = traj.position
        traj.correct_periodic = True
        pos_b = traj.position

        assert not np.allclose(pos_a, pos_b)
        assert snapshot_custom == pos_a
        traj.correct_periodic = False

    def test_position_at_frame(self, universe):
        traj = mn.entities.Trajectory(universe)
        traj.create_object()
        assert not np.allclose(traj._position_at_frame(1), traj._position_at_frame(3))

    @pytest.mark.parametrize(
        "correct,subframes,interpolate",
        itertools.product([True, False], [0, 1, 2, 3], [True, False]),
    )
    def test_mean_position(
        self, snapshot, subframes: int, correct: bool, interpolate: bool, universe
    ):
        traj = mn.entities.Trajectory(universe)
        traj.create_object()
        traj.correct_periodic = correct
        traj.subframes = subframes
        traj.interpolate = interpolate
        traj.subframes = 0
        assert np.allclose(traj.position_cache_mean(1), traj._position_at_frame(1))
        traj.average = 1
        assert not np.allclose(traj.position_cache_mean(1), traj._position_at_frame(1))
        assert snapshot == traj.position_cache_mean(1)
        assert snapshot == traj.cache

    def test_update_selection(self, snapshot_custom, universe):
        # to API add selections we currently have to operate on the UIList rather than the
        # universe itself, which isn't great

        traj = mn.entities.Trajectory(universe)
        traj.create_object()
        bpy.context.scene.frame_set(0)
        sel = traj.add_selection(
            name="custom_sel_1", selection_str="around 3.5 protein"
        )
        bpy.context.scene.frame_set(2)
        sel_1 = traj.named_attribute("custom_sel_1")
        bpy.context.scene.frame_set(4)
        sel_2 = traj.named_attribute("custom_sel_1")
        # when we are updating, the selection around the protein will change from frame
        # to frame
        assert not np.allclose(sel_1, sel_2)

        # if we stop the selection from updating, then even when we change the frame
        # the selection will remain the same
        sel.updating = False
        bpy.context.scene.frame_set(100)
        assert (sel_2 == traj.named_attribute("custom_sel_1")).all()
        # if we change the selection to updating, then the selection will be updated
        # and will no longer match with what came earlier
        sel.updating = False
        assert not (sel_2 != traj.named_attribute("custom_sel_1")).all()

    def test_selection_from_atomgroup(
        self,
        universe,
    ):
        ca_ag = universe.select_atoms("name CA")
        around_protein = universe.select_atoms("around 3.5 protein", updating=True)
        traj = mn.entities.Trajectory(universe)
        traj.create_object()
        bpy.context.scene.frame_set(0)
        traj.add_selection_from_atomgroup(atomgroup=ca_ag, name="ca")
        traj.add_selection_from_atomgroup(
            atomgroup=around_protein, name="around_protein"
        )
        bpy.context.scene.frame_set(2)
        sel_1 = traj.named_attribute("around_protein")
        bpy.context.scene.frame_set(4)
        sel_2 = traj.named_attribute("around_protein")
        assert not np.allclose(sel_1, sel_2)

    def test_save_persistance(
        self,
        snapshot_custom: NumpySnapshotExtension,
        tmp_path,
        universe,
        session: mn.session.MNSession,
    ):
        traj = mn.entities.Trajectory(universe)
        traj.create_object()
        traj.reset_playback()
        uuid = traj.uuid
        bpy.context.scene.frame_set(2)
        filepath = str(tmp_path / "test.blend")

        # test that we can save the file and it is created only after saving
        assert not os.path.exists(session.stashpath(filepath))
        bpy.ops.wm.save_as_mainfile(filepath=filepath)
        assert os.path.exists(filepath)
        assert os.path.exists(session.stashpath(filepath))
        del traj
        bpy.ops.wm.open_mainfile(filepath=filepath)

        traj = mn.session.get_session().trajectories[uuid]

        verts_frame_0 = traj.named_attribute("position")
        bpy.context.scene.frame_set(3)
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

    bpy.context.scene.frame_set(3)
    pos_b = traj.named_attribute("position")
    assert not np.allclose(pos_a, pos_b)

    for att in obj.data.attributes.keys():
        assert snapshot_custom == traj.named_attribute(att)
