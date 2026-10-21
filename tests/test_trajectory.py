import itertools
import os
import pickle
import tempfile
import bpy
import databpy
import MDAnalysis as mda
import numpy as np
import pytest
from numpy.testing import assert_almost_equal
import molecularnodes as mn
from .constants import data_dir
from .utils import NumpySnapshotExtension

pytestmark = [
    pytest.mark.filterwarnings("ignore:.*Empty string to select atoms.*:UserWarning"),
    pytest.mark.filterwarnings(
        "ignore:.*Unknown masses are set to.*:PendingDeprecationWarning"
    ),
]


def dummy_calculation_for_pickle_test(universe):
    """Module-level function for testing calculations pickling."""
    return universe.atoms.positions.mean(axis=0)


class TestMoleculeSingleStructure:
    @pytest.fixture
    def mol(self, fetch) -> mn.Molecule:
        return fetch("4ozs")

    def test_positions(self, mol: mn.Molecule):
        pos_0 = mol.position
        mol.position -= mol.centroid()
        pos_moved = mol.position
        assert not np.allclose(pos_0, pos_moved)
        # changing the scene frame calls set_frame() on all entities; a
        # single-frame structure should keep manually modified positions
        # rather than having them reset to the universe coordinates
        bpy.context.scene.frame_set(1)
        bpy.context.scene.frame_set(0)
        assert_almost_equal(pos_moved, mol.position)


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
        traj = mn.entities.Molecule(universe_with_bonds)
        assert len(traj.data.edges) == 28015

    def test_frames_to_collection(self, universe):
        traj = mn.entities.Molecule(universe, name="TestFramesToCollection")
        n_frames = traj.universe.trajectory.n_frames
        n_atoms = traj.universe.atoms.n_atoms

        # bake every frame
        frames = traj.frames_to_collection()
        assert len(frames.objects) == n_frames
        # each frame object has one vertex per atom
        for obj in frames.objects:
            assert len(obj.data.vertices) == n_atoms

        # a baked frame's positions match the universe at that frame
        traj.uframe = 2
        expected = traj._scaled_position
        got = databpy.named_attribute(
            frames.objects[f"{traj.name}_frame_0002"], "position"
        )
        assert np.allclose(got, expected)

        # a stepped range bakes fewer frames, and re-baking replaces (not appends)
        stepped = traj.frames_to_collection(start=0, stop=4, step=2)
        assert stepped is frames
        assert len(stepped.objects) == 2

        # step must be positive
        with pytest.raises(ValueError):
            traj.frames_to_collection(step=0)

    def test_attributes_added(self, universe):
        traj = mn.entities.Molecule(universe)
        attributes = traj.list_attributes()
        # check if all attributes are added.

        attribute_added = [
            "vdw_radii",
            "b_factor",
            "atomic_number",
            "res_id",
            "res_name",
            # "chain_id",
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
        traj = mn.entities.Molecule(universe, name="TestTrajectoryUpdate")
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
        traj = mn.entities.Molecule(universe)
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
        traj = mn.entities.Molecule(universe)
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
        traj = mn.entities.Molecule(univ_across_boundary)
        traj.subframes = 5
        bpy.context.scene.frame_set(2)
        pos_a = traj.position
        traj.correct_periodic = True
        pos_b = traj.position

        assert not np.allclose(pos_a, pos_b)
        assert snapshot_custom == pos_a
        traj.correct_periodic = False

    def test_correct_periodic_no_drift(self, univ_across_boundary):
        # the correction must only bridge the interpolation to the nearest
        # periodic image: on exact frame ticks the wrapped trajectory positions
        # are shown as-is. It used to mutate the cached positions, so each
        # frame was corrected against the previous already-shifted frame and
        # the trajectory unwrapped further out on every boundary crossing
        traj = mn.entities.Molecule(univ_across_boundary)
        traj.subframes = 3
        traj.interpolate = True
        traj.correct_periodic = True

        n_frames = traj.universe.trajectory.n_frames
        # walk the timeline through the subframes, as playback does
        for scene_frame in range((n_frames - 1) * 4 + 1):
            bpy.context.scene.frame_set(scene_frame)
            if scene_frame % 4 == 0:
                shown = traj.position.copy()
                raw = traj.frame_manager._position_at_frame(scene_frame // 4)
                assert np.allclose(shown, raw)
        traj.correct_periodic = False

    def test_correct_periodic_interpolates_averaged_frames(self, univ_across_boundary):
        # subframes between averaged frames must also bridge to the nearest
        # periodic image — the correction used to be skipped when averaging,
        # so a crossing atom swept across the whole box over the subframes
        traj = mn.entities.Molecule(univ_across_boundary)
        traj.subframes = 3
        traj.interpolate = True
        traj.average = 1
        traj.correct_periodic = True

        # within a segment each subframe step is a quarter of the (corrected)
        # frame-to-frame motion, far below the ~quarter-box sweeps of an
        # uncorrected crossing; steps onto exact ticks are skipped, as that is
        # where atoms legitimately wrap to the other side of the box
        dimensions = traj.universe.dimensions[:3] * traj.world_scale
        previous = None
        for scene_frame in range(4, 4 * 7 + 1):
            bpy.context.scene.frame_set(scene_frame)
            position = np.asarray(traj.position).copy()
            if previous is not None and scene_frame % 4 != 0:
                step = np.abs(position - previous)
                assert (step < dimensions * 0.2).all()
            previous = position
        traj.average = 0
        traj.correct_periodic = False

    def test_position_at_frame(self, universe):
        traj = mn.entities.Molecule(universe)
        assert not np.allclose(
            traj.frame_manager._position_at_frame(1),
            traj.frame_manager._position_at_frame(3),
        )

    @pytest.mark.parametrize(
        "correct,subframes,interpolate",
        list(itertools.product([True, False], [0, 1, 2, 3], [True, False])),
    )
    def test_mean_position(
        self, snapshot, subframes: int, correct: bool, interpolate: bool, universe
    ):
        traj = mn.entities.Molecule(universe)
        traj.correct_periodic = correct
        traj.subframes = subframes
        traj.interpolate = interpolate
        traj.subframes = 0
        traj.frame = 1
        pos_1 = traj.position
        assert np.allclose(traj.frame_manager.position_cache_mean(1), pos_1)
        traj.average = 1
        assert not np.allclose(
            traj.frame_manager.position_cache_mean(1),
            pos_1,
        )

    def test_manual_frame_raw(self, universe):
        # with no subframes or averaging, the manual frame is the universe frame
        traj = mn.entities.Molecule(universe)
        traj.update_with_scene = False
        traj.frame = 3
        assert np.allclose(traj.position, traj.frame_manager._position_at_frame(3))

    @pytest.mark.parametrize("interpolate", [True, False])
    def test_manual_frame_subframes(self, universe, interpolate: bool):
        # the manual frame steps through subframes like the scene frame does
        traj = mn.entities.Molecule(universe)
        traj.update_with_scene = False
        traj.subframes = 0
        traj.interpolate = interpolate
        traj.frame = 0
        verts_a = traj.position
        traj.frame = 1
        verts_b = traj.position
        assert not np.allclose(verts_a, verts_b)

        for subframes in [1, 2, 3, 4]:
            frame = 1
            fraction = frame % (subframes + 1) / (subframes + 1)
            traj.subframes = subframes
            traj.frame = frame
            verts_c = traj.position

            if interpolate:
                assert not np.allclose(verts_b, verts_c)
                assert np.allclose(verts_c, databpy.lerp(verts_a, verts_b, t=fraction))
            else:
                assert np.allclose(verts_a, verts_c)

    def test_manual_frame_average(self, universe):
        traj = mn.entities.Molecule(universe)
        traj.update_with_scene = False
        traj.frame = 1
        pos_1 = traj.position
        traj.average = 1
        expected = np.mean(
            [traj.frame_manager._position_at_frame(f) for f in [0, 1, 2]], axis=0
        )
        assert not np.allclose(traj.position, pos_1)
        assert np.allclose(traj.position, expected)

    def test_manual_frame_offset_ignored(self, universe):
        # the offset shifts where playback starts on the timeline, which has no
        # meaning for a manually chosen frame
        traj = mn.entities.Molecule(universe)
        traj.update_with_scene = False
        traj.frame = 2
        pos_2 = traj.position
        traj.offset = 1
        assert np.allclose(traj.position, pos_2)

    def test_manual_frame_clamp(self, universe):
        # the manual frame clamps to the last frame that maps onto the
        # trajectory, which is expanded by the subframes
        traj = mn.entities.Molecule(universe)
        traj.update_with_scene = False
        n_frames = universe.trajectory.n_frames
        traj.frame = 100
        assert traj.frame == n_frames - 1
        traj.subframes = 2
        traj.frame = 100
        assert traj.frame == (n_frames - 1) * 3
        assert np.allclose(
            traj.position, traj.frame_manager._position_at_frame(n_frames - 1)
        )

    def test_gui_selection_add_remove(self, universe):
        traj = mn.Molecule(universe)
        assert "selection_0" not in traj.list_attributes()
        bpy.ops.mn.trajectory_selection_add()
        assert "selection_0" in traj.list_attributes()
        assert traj.selections.ui_index == 0
        bpy.ops.mn.trajectory_selection_add()
        bpy.ops.mn.trajectory_selection_add()
        assert "selection_1" in traj.list_attributes()
        assert "selection_2" in traj.list_attributes()
        traj.selections.ui_index = 1
        bpy.ops.mn.trajectory_selection_remove()
        assert "selection_1" not in traj.list_attributes()
        assert traj.selections.ui_index == 1

        traj.selections.from_string("around 3.5 protein")
        traj.set_frame(0)
        sel_0 = traj["selection_1"]
        traj.set_frame(1)
        sel_1 = traj["selection_1"]
        assert not np.array_equal(sel_0, sel_1)
        traj.selections["selection_1"].updating = False
        traj.set_frame(2)
        assert np.array_equal(sel_1, traj["selection_1"])

        with pytest.raises(ValueError):
            traj.selections.remove("non_existent_selection")
        with pytest.raises(ValueError):
            traj.selections.remove(int)
        with pytest.raises(IndexError):
            traj.selections.remove(10)

    def test_update_selection(self, snapshot_custom, universe):
        # to API add selections we currently have to operate on the UIList rather than the
        # universe itself, which isn't great

        traj = mn.entities.Molecule(universe)
        bpy.context.scene.frame_set(0)
        sel = traj.selections.from_string("around 3.5 protein", name="custom_sel_1")
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
        traj = mn.entities.Molecule(universe)
        bpy.context.scene.frame_set(0)
        traj.selections.from_atomgroup(atomgroup=ca_ag, name="ca")
        traj.selections.from_atomgroup(atomgroup=around_protein, name="around_protein")
        bpy.context.scene.frame_set(2)
        sel_1 = traj.named_attribute("around_protein")
        bpy.context.scene.frame_set(4)
        sel_2 = traj.named_attribute("around_protein")
        assert not np.allclose(sel_1, sel_2)

    def test_save_persistance(
        self,
        snapshot_custom: NumpySnapshotExtension,
        tmp_path,
        session: mn.session.MNSession,
    ):
        # Load a new universe instead of a using the fixture as
        # this test modifies the filepaths in the universe
        top = data_dir / "md_ppr/box.gro"
        traj = data_dir / "md_ppr/first_5_frames.xtc"
        universe = mda.Universe(top, traj)
        traj = mn.entities.Molecule(universe)
        traj.reset_playback()
        uuid = traj.uuid
        bpy.context.scene.frame_set(2)

        # Add selections to test PyCapsule handling
        traj.selections.from_string("protein", name="protein_sel")
        traj.selections.from_string(
            "around 3.0 resid 1", name="dynamic_sel", updating=True
        )

        filepath = str(tmp_path / "test.blend")

        assert not os.path.exists(session.stashpath(filepath))
        try:
            bpy.ops.wm.save_as_mainfile(filepath=filepath)
        except Exception as e:
            if "cannot pickle 'PyCapsule' object" in str(e):
                pytest.fail(f"Session save failed with PyCapsule error: {e}")
            else:
                raise

        assert os.path.exists(filepath)
        assert os.path.exists(session.stashpath(filepath))
        # Rewind the trajectory to force a reopen of file by the reader
        # Fails if the trajectory paths aren't converted back to abs paths
        traj.universe.trajectory.rewind()
        del traj
        bpy.ops.wm.open_mainfile(filepath=filepath)

        traj = mn.session.get_session().get(uuid)

        # Verify selections were restored
        assert "protein_sel" in traj.list_attributes()
        assert "dynamic_sel" in traj.list_attributes()

        verts_frame_0 = traj.named_attribute("position")
        bpy.context.scene.frame_set(3)
        verts_frame_4 = traj.named_attribute("position")

        assert snapshot_custom == verts_frame_4
        assert not np.allclose(verts_frame_0, verts_frame_4)
        # Rewind the trajectory to force a reopen of file by the reader
        # Fails if the trajectory paths aren't abs paths after load
        traj.universe.trajectory.rewind()

    def test_get_view(self, universe):
        t1 = mn.Molecule(universe)
        # Note: When a frame is not specified, whatever is the current
        # universe frame should be used. In some random test failures
        # with bpy, the frame value was 3 for both cases leading to
        # an assertion failure. It is very likely a pytest issue.
        # Explicitly start with the frame reset to 0 for the universe
        universe.trajectory[0]
        # view of resid 1 at trajectory frame 0
        # view of resid 1
        v1 = t1.get_view(selection="resid 1")
        # view of resid 1 at trajectory frame 3
        v2 = t1.get_view(selection="resid 1", frame=3)
        assert v2 != v1
        # view of resid 2
        v3 = t1.get_view(selection="resid 2")
        assert v3 != v1 and v3 != v2
        # view of whole object
        v4 = t1.get_view()
        assert v4 != v1 and v4 != v2 and v4 != v3

    def test_gh_985(self):
        # ensure that we can create a Trajectory
        # from the SMILES representation of a simple molecule
        ethanol_smiles = "CCO"
        u = mda.Universe.from_smiles(ethanol_smiles)
        mn.Molecule(u)

    def test_trajectory_pickle_no_pycapsule_error(self, universe):
        """Test that Trajectory objects can be pickled without PyCapsule errors."""
        traj = mn.entities.Molecule(universe, name="TestPickleTrajectory")

        # Add selections and calculations
        traj.selections.from_string("protein", name="protein_sel")
        traj.selections.from_string(
            "around 5.0 resid 1", name="around_sel", updating=True
        )
        traj.calculations["center_of_mass"] = dummy_calculation_for_pickle_test
        # Change frame to check frame restore
        traj.uframe = 1

        with tempfile.NamedTemporaryFile() as tmp_file:
            try:
                pickle.dump(traj, tmp_file)
                tmp_file.flush()

                tmp_file.seek(0)
                restored_traj = pickle.load(tmp_file)

                # Verify restoration
                assert restored_traj.name == "TestPickleTrajectory"
                assert hasattr(restored_traj, "universe")
                assert hasattr(restored_traj, "frame_manager")
                assert hasattr(restored_traj, "selections")
                assert hasattr(restored_traj, "annotations")
                assert hasattr(restored_traj, "calculations")

                # Verify universe restoration
                assert hasattr(restored_traj.universe, "atoms")
                assert restored_traj.universe.atoms.n_atoms == universe.atoms.n_atoms
                assert restored_traj.universe.trajectory.frame == 1

                # Verify circular references
                assert restored_traj.frame_manager.trajectory is restored_traj
                assert restored_traj.selections.trajectory is restored_traj

            except TypeError as e:
                if "cannot pickle 'PyCapsule' object" in str(e):
                    pytest.fail(f"PyCapsule pickle error not fixed: {e}")
                else:
                    raise

    def test_trajectory_pickle_deserialization_failure(self, universe):
        """Test that trajectory deserialization fails fast with clear errors when files are missing."""
        traj = mn.entities.Molecule(universe, name="TestFailedDeserialization")

        # Pickle the trajectory
        with tempfile.NamedTemporaryFile() as tmp_file:
            pickle.dump(traj, tmp_file)
            tmp_file.flush()
            tmp_file.seek(0)

            # Manually corrupt the pickled data to simulate missing/invalid file paths
            # by modifying the pickled trajectory to reference non-existent files
            pickled_data = pickle.load(tmp_file)

            # Create a new state dict with invalid file paths
            state = pickled_data.__getstate__()
            state["_universe_topology"] = "/nonexistent/path/topology.pdb"
            state["_universe_trajectory"] = "/nonexistent/path/trajectory.xtc"

            # Attempt to restore - this should raise RuntimeError, not create broken object
            new_traj = mn.entities.Molecule.__new__(mn.entities.Molecule)
            with pytest.raises(RuntimeError) as exc_info:
                new_traj.__setstate__(state)

            # Verify the error message is descriptive
            error_msg = str(exc_info.value)
            assert "Failed to restore Molecule from saved session" in error_msg
            assert "Could not recreate MDAnalysis Universe" in error_msg
            assert "/nonexistent/path/topology.pdb" in error_msg
            assert "/nonexistent/path/trajectory.xtc" in error_msg

    def test_dssp(self, snapshot, universe):
        initial_frame = universe.trajectory.frame
        t = mn.Molecule(universe).add_style("cartoon")
        # sec_struct is computed once from the current frame on import, so the
        # cartoon has structure to show before the dssp module is initialized
        initial_sec_struct = t["sec_struct"]
        assert np.any(initial_sec_struct == 1)
        assert "sec_struct" not in t.calculations
        # initialize dssp
        t.dssp.init()
        # per-frame dssp on the import frame matches the import-time result
        t.set_frame(initial_frame)
        assert np.array_equal(t["sec_struct"], initial_sec_struct)
        # default dssp is per-frame
        t.set_frame(1)
        frame_sec_struct = t["sec_struct"]
        # change dssp to none
        t.dssp.show_none()
        no_sec_struct = t["sec_struct"]
        assert not np.allclose(frame_sec_struct, no_sec_struct)
        assert snapshot == no_sec_struct
        # change dssp to sliding-window-average
        t.dssp.show_sliding_window_average(window_size=5)
        sw_sec_struct = t["sec_struct"]
        assert snapshot == sw_sec_struct
        # change dssp to trajectory-average
        t.dssp.show_trajectory_average(threshold=0.5)
        avg_sec_struct = t["sec_struct"]
        assert snapshot == avg_sec_struct
        # per-frame, sliding-window-average and trajectory-average
        # should all be different from each other
        assert not np.allclose(frame_sec_struct, sw_sec_struct)
        assert not np.allclose(sw_sec_struct, avg_sec_struct)
        # per-frame, sliding-window-average and trajectory-average
        # should all be different from none
        assert not np.allclose(no_sec_struct, frame_sec_struct)
        assert not np.allclose(no_sec_struct, sw_sec_struct)
        assert not np.allclose(no_sec_struct, avg_sec_struct)

    def test_dssp_selection(self, universe):
        t = mn.Molecule(universe).add_style("cartoon")
        t.dssp.init(selection="protein and resid 1-100")
        t.set_frame(1)
        sec_struct = t["sec_struct"]
        resids = universe.atoms.resids
        is_protein = np.isin(
            universe.atoms.indices, universe.select_atoms("protein").indices
        )
        selected = is_protein & (resids <= 100)
        # selected residues are assigned helix / sheet / loop
        assert np.all(np.isin(sec_struct[selected], [1, 2, 3]))
        assert np.any(sec_struct[selected] < 3)
        # protein residues outside the selection fall back to loop
        assert np.all(sec_struct[is_protein & ~selected] == 3)
        # non-protein atoms get no secondary structure
        assert np.all(sec_struct[~is_protein] == 0)
        # a default run assigns structure beyond the selection
        t_full = mn.Molecule(universe).add_style("cartoon")
        t_full.dssp.init()
        t_full.set_frame(1)
        assert np.any(t_full["sec_struct"][is_protein & ~selected] < 3)

    def test_dssp_set_selection(self, universe):
        t = mn.Molecule(universe).add_style("cartoon")
        t.dssp.init(selection="protein and resid 1-100")
        assert t.props.dssp.selection == "protein and resid 1-100"
        t.set_frame(1)
        resids = universe.atoms.resids
        is_protein = np.isin(
            universe.atoms.indices, universe.select_atoms("protein").indices
        )
        excluded = is_protein & (resids > 100)
        assert np.all(t["sec_struct"][excluded] == 3)
        # widening the selection part way through assigns structure to the rest
        t.dssp.set_selection("protein")
        assert t.props.dssp.selection == "protein"
        t.set_frame(1)
        assert np.any(t["sec_struct"][excluded] < 3)
        # a bad selection raises and leaves the previous selection working
        with pytest.raises(Exception):
            t.dssp.set_selection("protein and bogus")
        t.set_frame(2)
        assert np.any(t["sec_struct"][excluded] < 3)
        # the display option survives changing the selection
        t.dssp.show_trajectory_average(threshold=0.5)
        t.dssp.set_selection("protein and resid 1-100")
        assert t.dssp._display_option == "trajectory-average"
        t.set_frame(1)
        assert np.all(t["sec_struct"][excluded] == 3)
        assert np.any(t["sec_struct"][is_protein & ~excluded] < 3)

    def test_reload_operator(self, universe):
        traj = mn.entities.Molecule(universe)
        traj.add_style("cartoon")
        obj_name = traj.name
        session = mn.session.get_session()
        session.entities.pop(traj.uuid)

        bpy.data.objects[obj_name].select_set(True)
        bpy.ops.mn.session_reload_item()

        assert len(session.entities) == 1

    @staticmethod
    def _universe_in_memory(velocities: bool, forces: bool) -> mda.Universe:
        # the test trajectory has neither, so give a copy of it random
        # per-frame velocities and forces in memory
        u = mda.Universe(
            data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
        )
        coords = np.array([ts.positions.copy() for ts in u.trajectory])
        rng = np.random.default_rng(0)
        kwargs = {}
        if velocities:
            kwargs["velocities"] = rng.normal(size=coords.shape).astype(np.float32)
        if forces:
            kwargs["forces"] = rng.normal(size=coords.shape).astype(np.float32) * 100
        u.load_new(coords, format=mda.coordinates.memory.MemoryReader, **kwargs)
        return u

    @pytest.mark.parametrize(
        "velocities,forces", [(True, True), (True, False), (False, True)]
    )
    def test_velocity_force_attributes(self, velocities: bool, forces: bool):
        traj = mn.entities.Molecule(self._universe_in_memory(velocities, forces))
        names = traj.list_attributes()
        assert ("velocity" in names) == velocities
        assert ("force" in names) == forces

        for frame in [0, 2, 4]:
            bpy.context.scene.frame_set(frame)
            ts = traj.universe.trajectory[frame]
            if velocities:
                assert np.allclose(
                    traj.named_attribute("velocity"), ts.velocities * traj.world_scale
                )
            if forces:
                assert np.allclose(traj.named_attribute("force"), ts.forces)

    def test_velocity_force_absent(self, universe):
        traj = mn.entities.Molecule(universe)
        bpy.context.scene.frame_set(1)
        names = traj.list_attributes()
        assert "velocity" not in names
        assert "force" not in names

    def test_velocity_force_subframes(self):
        traj = mn.entities.Molecule(self._universe_in_memory(True, True))
        traj.subframes = 3
        traj.interpolate = True
        u = traj.universe
        vel_a = u.trajectory[1].velocities * traj.world_scale
        force_a = u.trajectory[1].forces.copy()
        vel_b = u.trajectory[2].velocities * traj.world_scale
        force_b = u.trajectory[2].forces.copy()

        # with 3 subframes, scene frame 6 is halfway from frame 1 to frame 2
        bpy.context.scene.frame_set(6)
        assert np.allclose(
            traj.named_attribute("velocity"), databpy.lerp(vel_a, vel_b, t=0.5)
        )
        assert np.allclose(
            traj.named_attribute("force"), databpy.lerp(force_a, force_b, t=0.5)
        )

        # without interpolation subframes hold the previous frame's values
        traj.interpolate = False
        bpy.context.scene.frame_set(7)
        assert np.allclose(traj.named_attribute("velocity"), vel_a)

    def test_velocity_force_average(self):
        traj = mn.entities.Molecule(self._universe_in_memory(True, True))
        traj.average = 1
        u = traj.universe
        expected = np.mean([u.trajectory[f].velocities for f in (1, 2, 3)], axis=0)
        expected_force = np.mean([u.trajectory[f].forces for f in (1, 2, 3)], axis=0)
        bpy.context.scene.frame_set(2)
        assert np.allclose(
            traj.named_attribute("velocity"), expected * traj.world_scale
        )
        assert np.allclose(traj.named_attribute("force"), expected_force)

    def _write_trr(self, path, velocity_frames=(), force_frames=()) -> mda.Universe:
        # GROMACS writes velocities and forces at their own strides, so a TRR
        # file can have them on only some frames
        source = self._universe_in_memory(True, True)
        with mda.Writer(str(path), n_atoms=source.atoms.n_atoms) as writer:
            for ts in source.trajectory:
                ts.has_velocities = ts.frame in velocity_frames
                ts.has_forces = ts.frame in force_frames
                writer.write(source)
        return mda.Universe(data_dir / "md_ppr/box.gro", str(path))

    def test_velocity_force_trr_gaps(self, tmp_path):
        # velocities on frames 1 and 3, forces on frame 2, neither on frame 0
        u = self._write_trr(tmp_path / "gaps.trr", (1, 3), (2,))
        reference = mda.Universe(
            data_dir / "md_ppr/box.gro", str(tmp_path / "gaps.trr")
        )
        vel = {f: reference.trajectory[f].velocities.copy() for f in (1, 3)}
        force = reference.trajectory[2].forces.copy()

        traj = mn.entities.Molecule(u)
        bpy.context.scene.frame_set(0)
        assert "velocity" not in traj.list_attributes()
        assert "force" not in traj.list_attributes()
        assert traj.object.mn.frame_attributes == "position"

        # an attribute appears on the first frame that has it, and frames
        # without it keep the values last stored
        shown_vel = {1: 1, 2: 1, 3: 3, 4: 3, 0: 3}
        for frame in [1, 2, 3, 4, 0]:
            bpy.context.scene.frame_set(frame)
            assert traj.uframe == frame
            assert np.allclose(
                traj.named_attribute("velocity"),
                vel[shown_vel[frame]] * traj.world_scale,
            )
            if frame == 1:
                assert "force" not in traj.list_attributes()
            else:
                assert np.allclose(traj.named_attribute("force"), force)
        assert traj.object.mn.frame_attributes == "position,velocity,force"

        # interpolating towards a frame without velocities uses this frame's:
        # scene frame 2 is halfway from frame 1 to frame 2
        traj.subframes = 1
        traj.interpolate = True
        bpy.context.scene.frame_set(2)
        assert np.allclose(traj.named_attribute("velocity"), vel[1] * traj.world_scale)
        assert np.allclose(traj.named_attribute("force"), force)

    def test_frame_attribute_info(self, tmp_path):
        # velocities on frames 1 and 3, forces on frame 2
        traj = mn.entities.Molecule(
            self._write_trr(tmp_path / "info.trr", (1, 3), (2,))
        )
        for frame in [1, 2]:
            bpy.context.scene.frame_set(frame)

        info = traj.frame_attribute_info("velocity")
        assert info["source"] == "ts.velocities"
        assert info["units"] == "Blender units/ps (Å/ps × 0.1)"
        # frame 2 has no velocities, so frame 1's are held
        assert info["frame"] == 1
        assert not info["current"]
        lengths = np.linalg.norm(traj.named_attribute("velocity"), axis=1)
        assert np.allclose(info["range"], (lengths.min(), lengths.max()))

        info = traj.frame_attribute_info("force")
        assert info["source"] == "ts.forces"
        assert info["units"] == "kJ/(mol·Å)"
        assert info["frame"] == 2
        assert info["current"]

        assert traj.frame_attribute_info("position")["source"] == "ts.positions"

        # the range is recomputed once the attributes are next stored
        bpy.context.scene.frame_set(3)
        info = traj.frame_attribute_info("velocity")
        lengths = np.linalg.norm(traj.named_attribute("velocity"), axis=1)
        assert info["current"]
        assert np.allclose(info["range"], (lengths.min(), lengths.max()))

    def test_velocity_force_trr_positions_only(self, tmp_path):
        traj = mn.entities.Molecule(self._write_trr(tmp_path / "x.trr"))
        bpy.context.scene.frame_set(1)
        assert "velocity" not in traj.list_attributes()
        assert "force" not in traj.list_attributes()

    def test_ts_data_pdb_models(self, tmp_path):
        # each model of a multi-model PDB file has its own B-factors and
        # occupancies, which MDAnalysis reads into ts.data
        source = mda.Universe(
            data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
        )
        atoms = source.atoms[:50]
        source.add_TopologyAttr("tempfactors")
        source.add_TopologyAttr("occupancies")
        path = tmp_path / "models.pdb"
        rng = np.random.default_rng(0)
        b_factors = rng.uniform(0, 100, size=(3, 50)).round(2)
        occupancies = rng.uniform(0, 1, size=(3, 50)).round(2)
        with mda.Writer(str(path), multiframe=True) as writer:
            for i, _ in enumerate(source.trajectory[:3]):
                atoms.tempfactors = b_factors[i]
                atoms.occupancies = occupancies[i]
                writer.write(atoms)

        traj = mn.entities.Molecule(mda.Universe(str(path)))
        assert traj.object.mn.frame_attributes == "position,occupancy,b_factor"
        for frame in range(3):
            bpy.context.scene.frame_set(frame)
            assert np.allclose(traj.named_attribute("b_factor"), b_factors[frame])
            assert np.allclose(traj.named_attribute("occupancy"), occupancies[frame])

        traj.subframes = 1
        traj.interpolate = True
        bpy.context.scene.frame_set(1)
        assert np.allclose(
            traj.named_attribute("b_factor"), b_factors[:2].mean(axis=0), atol=1e-4
        )
        assert np.allclose(
            traj.named_attribute("occupancy"),
            occupancies[:2].mean(axis=0),
            atol=1e-4,
        )

    @staticmethod
    def _write_models_pdb(path, b_factors, occupancies):
        # a multi-model PDB file of 50 atoms with these values for each model
        source = mda.Universe(
            data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
        )
        atoms = source.atoms[:50]
        source.add_TopologyAttr("tempfactors")
        source.add_TopologyAttr("occupancies")
        with mda.Writer(str(path), multiframe=True) as writer:
            for i, _ in enumerate(source.trajectory[: len(b_factors)]):
                atoms.tempfactors = b_factors[i]
                atoms.occupancies = occupancies[i]
                writer.write(atoms)

    @staticmethod
    def _assert_models_follow_frames(mol, b_factors, occupancies):
        for frame in range(len(b_factors)):
            bpy.context.scene.frame_set(frame)
            assert np.allclose(mol.named_attribute("b_factor"), b_factors[frame])
            assert np.allclose(mol.named_attribute("occupancy"), occupancies[frame])

    def test_ts_data_pdb_models_from_file(self, tmp_path):
        # a multi-model PDB file loaded on its own is parsed by biotite, whose
        # stack has the first model's values, and each frame carries its own
        # model's B-factors and occupancies in ts.data
        rng = np.random.default_rng(0)
        b_factors = rng.uniform(0, 100, size=(3, 50)).round(2)
        occupancies = rng.uniform(0, 1, size=(3, 50)).round(2)
        path = tmp_path / "models.pdb"
        self._write_models_pdb(path, b_factors, occupancies)

        mol = mn.Molecule.load(path)
        assert mol.object.mn.frame_attributes == "position,occupancy,b_factor"
        self._assert_models_follow_frames(mol, b_factors, occupancies)

        mol.subframes = 1
        mol.interpolate = True
        bpy.context.scene.frame_set(1)
        assert np.allclose(
            mol.named_attribute("b_factor"), b_factors[:2].mean(axis=0), atol=1e-4
        )
        assert np.allclose(
            mol.named_attribute("occupancy"), occupancies[:2].mean(axis=0), atol=1e-4
        )

        mol.subframes = 0
        mol.average = 1
        bpy.context.scene.frame_set(1)
        assert np.allclose(
            mol.named_attribute("b_factor"), b_factors.mean(axis=0), atol=1e-4
        )

    def test_ts_data_pdb_models_from_file_identical(self, tmp_path):
        # values that are the same in every model are not per-frame attributes
        rng = np.random.default_rng(0)
        b_factors = np.tile(rng.uniform(0, 100, size=50).round(2), (3, 1))
        occupancies = rng.uniform(0, 1, size=(3, 50)).round(2)
        path = tmp_path / "models.pdb"
        self._write_models_pdb(path, b_factors, occupancies)

        mol = mn.Molecule.load(path)
        assert "tempfactor" not in mol.universe.trajectory.ts.data
        assert mol.object.mn.frame_attributes == "position,occupancy"
        self._assert_models_follow_frames(mol, b_factors, occupancies)

    def test_ts_data_cif_models_from_file(self, tmp_path):
        import biotite.structure.io.pdb as pdb
        import biotite.structure.io.pdbx as pdbx

        rng = np.random.default_rng(1)
        b_factors = rng.uniform(0, 100, size=(3, 50)).round(2)
        occupancies = rng.uniform(0, 1, size=(3, 50)).round(2)
        pdb_path = tmp_path / "models.pdb"
        self._write_models_pdb(pdb_path, b_factors, occupancies)
        stack = pdb.get_structure(pdb.PDBFile.read(pdb_path))
        cif = pdbx.CIFFile()
        pdbx.set_structure(cif, stack)
        # the rows of `atom_site` are the atoms of each model in turn
        atom_site = cif.block["atom_site"]
        atom_site["B_iso_or_equiv"] = b_factors.ravel()
        atom_site["occupancy"] = occupancies.ravel()
        path = tmp_path / "models.cif"
        cif.write(path)

        mol = mn.Molecule.load(path)
        assert mol.object.mn.frame_attributes == "position,occupancy,b_factor"
        self._assert_models_follow_frames(mol, b_factors, occupancies)

    def test_ts_data_pdb_models_from_file_session(self, tmp_path, session):
        from molecularnodes.entities.reload import reload_entity

        rng = np.random.default_rng(2)
        b_factors = rng.uniform(0, 100, size=(3, 50)).round(2)
        occupancies = rng.uniform(0, 1, size=(3, 50)).round(2)
        path = tmp_path / "models.pdb"
        self._write_models_pdb(path, b_factors, occupancies)

        mol = mn.Molecule.load(path)
        uuid = mol.uuid
        bpy.context.scene.frame_set(2)
        filepath = str(tmp_path / "models.blend")
        bpy.ops.wm.save_as_mainfile(filepath=filepath)
        del mol
        bpy.ops.wm.open_mainfile(filepath=filepath)

        # the universe is rebuilt from the file, with each model's values
        mol = mn.session.get_session().get(uuid)
        assert "tempfactor" in mol.universe.trajectory.ts.data
        self._assert_models_follow_frames(mol, b_factors, occupancies)

        # as is an entity reloaded from its object
        obj = mol.object
        mn.session.get_session().remove_entity(uuid)
        mol = reload_entity(obj)
        self._assert_models_follow_frames(mol, b_factors, occupancies)

    def test_ts_data_lammps_dump(self, tmp_path):
        # extra columns of a LAMMPS dump are read into ts.data, and `q` updates
        # the charges
        charges = np.array([[0.1, -0.2, 0.3, -0.4], [0.5, -0.6, 0.7, -0.8]])
        lines = []
        for step, q in enumerate(charges):
            lines += [
                "ITEM: TIMESTEP",
                str(step),
                "ITEM: NUMBER OF ATOMS",
                "4",
                "ITEM: BOX BOUNDS pp pp pp",
                "0 10",
                "0 10",
                "0 10",
                "ITEM: ATOMS id type x y z q",
            ]
            lines += [f"{i + 1} 1 {i + step} 1.0 1.0 {q[i]}" for i in range(4)]
        path = tmp_path / "charges.lammpsdump"
        path.write_text("\n".join(lines) + "\n")

        u = mda.Universe(str(path), format="LAMMPSDUMP", additional_columns=["q"])
        traj = mn.entities.Molecule(u)
        for frame in range(2):
            bpy.context.scene.frame_set(frame)
            assert np.allclose(traj.named_attribute("charge"), charges[frame])
        assert "q" not in traj.list_attributes()

    def test_ts_data_integer_not_interpolated(self):
        # per-atom integer values hold the current frame's values on subframes
        # a file reader, as the MemoryReader applies transformations once up front
        u = mda.Universe(
            data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
        )
        n_atoms = u.atoms.n_atoms

        def add_flag(ts):
            ts.data["flag"] = np.full(n_atoms, ts.frame * 10, dtype=np.int32)
            return ts

        u.trajectory.add_transformations(add_flag)
        traj = mn.entities.Molecule(u)
        traj.subframes = 1
        traj.interpolate = True
        # scene frame 3 is halfway from frame 1 (10) to frame 2 (20)
        bpy.context.scene.frame_set(3)
        assert np.all(traj.named_attribute("flag") == 10)
        bpy.context.scene.frame_set(4)
        assert np.all(traj.named_attribute("flag") == 20)

    def test_frame_attributes_xtc(self, universe):
        traj = mn.entities.Molecule(universe)
        bpy.context.scene.frame_set(1)
        assert traj.object.mn.frame_attributes == "position"

    def test_frame_attributes_single_frame(self, universe_with_bonds):
        traj = mn.entities.Molecule(universe_with_bonds)
        assert traj.object.mn.frame_attributes == ""

    def test_ts_data_builtin_collision(self, tmp_path):
        # a LAMMPS dump column named like a built-in attribute is stored with a
        # `ts_` prefix and leaves the built-in attribute alone
        lines = []
        for step in range(2):
            lines += [
                "ITEM: TIMESTEP",
                str(step),
                "ITEM: NUMBER OF ATOMS",
                "4",
                "ITEM: BOX BOUNDS pp pp pp",
                "0 10",
                "0 10",
                "0 10",
                "ITEM: ATOMS id type x y z mass",
            ]
            lines += [f"{i + 1} 1 {i + step} 1.0 1.0 {100 + step}" for i in range(4)]
        path = tmp_path / "mass.lammpsdump"
        path.write_text("\n".join(lines) + "\n")

        u = mda.Universe(str(path), format="LAMMPSDUMP", additional_columns=["mass"])
        traj = mn.entities.Molecule(u)
        mass = traj.named_attribute("mass").copy()
        for frame in range(2):
            bpy.context.scene.frame_set(frame)
            assert np.allclose(traj.named_attribute("ts_mass"), 100 + frame)
            assert np.array_equal(traj.named_attribute("mass"), mass)
        assert traj.object.mn.frame_attributes == "position,ts_mass"
        assert traj.frame_attribute_info("ts_mass")["source"] == 'ts.data["mass"]'
        assert traj.frame_attribute_info("ts_mass")["units"] == ""


@pytest.mark.parametrize("topology", ["pent/prot_ion.tpr", "pent/TOPOL2.pdb"])
def test_martini(snapshot, topology):
    universe = mda.Universe(
        data_dir / "martini" / topology, data_dir / "martini/pent/PENT2_100frames.xtc"
    )
    traj = mn.entities.Molecule(universe)
    bpy.context.scene.frame_set(0)
    pos_a = traj.named_attribute("position")

    bpy.context.scene.frame_set(3)
    pos_b = traj.named_attribute("position")
    assert not np.allclose(pos_a, pos_b)

    for att in traj.list_attributes():
        assert snapshot == traj[att]
