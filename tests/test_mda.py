import bpy
import os
import pytest
import molecularnodes as mn
from . import utils
from molecularnodes.io.md import HAS_mda
from molecularnodes.blender.obj import get_attribute

if HAS_mda:
    import MDAnalysis as mda
import numpy as np
from .constants import data_dir
from .utils import remove_all_molecule_objects, sample_attribute, NumpySnapshotExtension


@pytest.mark.skipif(not HAS_mda, reason="MDAnalysis is not installed")
class TestMDA:
    @pytest.fixture(scope="module")
    def mda_session(self):
        mda_session = mn.io.MDAnalysisSession()
        return mda_session

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

    def test_persistent_handlers_added(self, mda_session):
        assert bpy.app.handlers.load_post[-1].__name__ == "_rejuvenate_universe"
        assert bpy.app.handlers.save_post[-1].__name__ == "_sync_universe"

    def test_create_mda_session(self, mda_session):
        assert mda_session is not None
        assert mda_session.world_scale == 0.01

    def reload_mda_session(self, mda_session):
        with pytest.warns(UserWarning, match="The existing mda session"):
            mda_session_2 = mn.mda.create_session()

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_show_universe(
        self, snapshot_custom: NumpySnapshotExtension, in_memory, mda_session, universe
    ):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe, in_memory=in_memory)
        bob = bpy.data.objects["atoms"]
        assert snapshot_custom == sample_attribute(bob, "position")

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_same_name_atoms(self, in_memory, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe, in_memory=in_memory)

        with pytest.warns(UserWarning, match="The name of the object is changed"):
            mda_session.show(universe, in_memory=in_memory)

        bob_1 = bpy.data.objects["atoms"]
        bob_2 = bpy.data.objects["atoms.001"]
        verts_1 = mn.blender.obj.get_attribute(bob_1, "position")
        verts_2 = mn.blender.obj.get_attribute(bob_2, "position")

        assert (verts_1 == verts_2).all()

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_show_multiple_selection(
        self, snapshot_custom: NumpySnapshotExtension, in_memory, mda_session, universe
    ):
        remove_all_molecule_objects(mda_session)
        custom_selections = {"name_ca": "name CA"}
        mda_session.show(
            universe,
            in_memory=in_memory,
            name="protein",
            selection="protein",
            custom_selections=custom_selections,
        )
        bob = bpy.data.objects["protein"]
        assert snapshot_custom == sample_attribute(bob, "posiiton")

        # different bahavior in_memory or not.
        if not in_memory:
            bob_ca = bpy.data.objects["name_ca"]
            assert snapshot_custom == sample_attribute(bob_ca, "position")
        else:
            # attribute is added as name_ca.
            assert "name_ca" in bob.data.attributes.keys()

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_include_bonds(self, in_memory, mda_session, universe_with_bonds):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe_with_bonds, in_memory=in_memory)
        obj = bpy.data.objects["atoms"]
        assert obj.data.edges.items() != []

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_attributes_added(self, in_memory, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe, in_memory=in_memory)
        obj = bpy.data.objects["atoms"]
        attributes = obj.data.attributes.keys()
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

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_trajectory_update(
        self, snapshot_custom: NumpySnapshotExtension, in_memory, universe
    ):
        # remove_all_molecule_objects(mda_session)
        mda_session = mn.io.MDAnalysisSession()

        obj = mda_session.show(universe, in_memory=in_memory, style="ribbon")
        node = mn.blender.nodes.get_style_node(obj)
        group = obj.modifiers["MolecularNodes"].node_group
        if in_memory:
            node = group.nodes["MN_animate_value"]
            node.inputs["Frame: Start"].default_value = 0
            node.inputs["Frame: End"].default_value = 4

        if "EEVEE" in node.inputs.keys():
            node.inputs["EEVEE"].default_value = True

        if in_memory:
            mn.blender.nodes.realize_instances(obj)

        n = 100

        pos_a = sample_attribute(obj, "position", n=n, evaluate=in_memory)
        assert snapshot_custom == pos_a

        # change blender frame to 4
        next_frame = 200 if in_memory else 4
        bpy.context.scene.frame_set(next_frame)

        # if in_memory:
        #     socket.default_value = 250
        pos_b = sample_attribute(obj, "position", n=n, evaluate=in_memory)
        assert snapshot_custom == pos_b

        assert not np.isclose(pos_a, pos_b).all()

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_show_updated_atoms(
        self, snapshot_custom: NumpySnapshotExtension, in_memory, mda_session, universe
    ):
        remove_all_molecule_objects(mda_session)
        updating_ag = universe.select_atoms("around 5 resid 1", updating=True)
        mda_session.show(updating_ag, in_memory=in_memory, style="vdw")

        bob = bpy.data.objects["atoms"]
        nodes = bob.modifiers["MolecularNodes"].node_group.nodes
        for node in nodes:
            for input in node.inputs:
                if input.name == "Frame: Start":
                    input.default_value = 0
                elif input.name == "Frame: End":
                    input.default_value = 4
                elif input.name == "To Max":
                    input.default_value = 4
                elif input.name == "EEVEE":
                    input.default_value = True

        mn.blender.nodes.realize_instances(bob)

        bpy.context.scene.frame_set(0)
        verts_frame_0 = get_attribute(bob, "position", evaluate=True)
        assert snapshot_custom == verts_frame_0

        # change blender frame to 1
        bpy.context.scene.frame_set(1)
        # bob = bpy.data.objects["atoms"]
        verts_frame_1 = get_attribute(bob, "position", evaluate=True)

        assert snapshot_custom == verts_frame_1

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_update_deleted_objects(self, in_memory, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe, in_memory=in_memory)
        bpy.data.objects.remove(bpy.data.objects["atoms"])

        # trigger depsgraph_update_post handler
        # by creating a new object
        bpy.ops.mesh.primitive_cube_add()

        assert mda_session.universe_reps == {}
        assert mda_session.atom_reps == {}
        assert mda_session.rep_names == []

        # remove the cube
        bpy.data.objects.remove(bpy.data.objects["Cube"])

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_save_persistance(
        self,
        snapshot_custom: NumpySnapshotExtension,
        tmp_path,
        in_memory,
        mda_session,
        universe,
    ):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe, in_memory=in_memory)
        # save
        bpy.ops.wm.save_as_mainfile(filepath=str(tmp_path / "test.blend"))

        assert os.path.exists(str(tmp_path / "test.mda_session"))

        # reload
        remove_all_molecule_objects(mda_session)
        bpy.ops.wm.open_mainfile(filepath=str(tmp_path / "test.blend"))
        bob = bpy.data.objects["atoms"]
        verts_frame_0 = mn.blender.obj.get_attribute(bob, "position")

        # change blender frame to 1
        bpy.context.scene.frame_set(1)
        bob = bpy.data.objects["atoms"]
        verts_frame_1 = mn.blender.obj.get_attribute(bob, "position")
        assert snapshot_custom == verts_frame_1

        assert not np.isclose(verts_frame_0, verts_frame_1).all()


@pytest.mark.skipif(not HAS_mda, reason="MDAnalysis is not installed")
class TestMDA_FrameMapping:
    @pytest.fixture(scope="module")
    def mda_session(self):
        mda_session = mn.io.MDAnalysisSession()
        return mda_session

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

    def test_persistent_handlers_added(self, mda_session):
        assert bpy.app.handlers.load_post[-1].__name__ == "_rejuvenate_universe"
        assert bpy.app.handlers.save_post[-1].__name__ == "_sync_universe"

    def test_create_mda_session(self, mda_session):
        assert mda_session is not None
        assert mda_session.world_scale == 0.01

    def reload_mda_session(self, mda_session):
        with pytest.warns(UserWarning, match="The existing mda session"):
            mda_session_2 = mn.mda.create_session()

    def test_frame_mapping(self, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        obj = mda_session.show(universe, frame_mapping=[0, 0, 1, 2, 4])

        bpy.context.scene.frame_set(0)
        verts_a = get_attribute(obj, "position")

        bpy.context.scene.frame_set(1)
        verts_b = get_attribute(obj, "position")
        # test the frame mapping works, that nothing has changed becuase of the mapping
        assert np.isclose(verts_a, verts_b).all()

        bpy.context.scene.frame_set(2)
        verts_b = get_attribute(obj, "position")
        # test that something has now changed
        assert not np.isclose(verts_a, verts_b).all()

    def test_subframes(self, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe)

        obj = bpy.data.objects["atoms"]
        bpy.context.scene.frame_set(0)
        verts_a = get_attribute(obj, "position")

        bpy.context.scene.frame_set(1)
        verts_b = get_attribute(obj, "position")
        # should be no difference because not using subframes
        assert not np.isclose(verts_a, verts_b).all()

        for subframes in [1, 2, 3, 4]:
            frame = 1
            fraction = frame % (subframes + 1) / (subframes + 1)
            obj.mn["subframes"] = subframes
            bpy.context.scene.frame_set(frame)
            verts_c = get_attribute(obj, "position")
            # now using subframes, there should be a difference
            assert not np.isclose(verts_b, verts_c).all()

            assert np.isclose(
                verts_c, mn.utils.lerp(verts_a, verts_b, t=fraction)
            ).all()

    def test_subframe_mapping(self, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe, in_memory=False, frame_mapping=[0, 0, 1, 2, 3])

        bob = bpy.data.objects["atoms"]
        bpy.context.scene.frame_set(0)
        verts_a = get_attribute(bob, "position")

        bpy.context.scene.frame_set(1)
        verts_b = get_attribute(bob, "position")
        assert np.isclose(verts_a, verts_b).all()

        bpy.context.scene.frame_set(2)
        verts_b = get_attribute(bob, "position")
        assert not np.isclose(verts_a, verts_b).all()

        bob.mn["subframes"] = 1
        bpy.context.scene.frame_set(3)
        verts_c = get_attribute(bob, "position")

        assert not np.isclose(verts_b, verts_c).all()
        assert np.isclose(verts_c, mn.utils.lerp(verts_a, verts_b, 0.5)).all()


@pytest.mark.parametrize("toplogy", ["pent/prot_ion.tpr", "pent/TOPOL2.pdb"])
def test_martini(snapshot_custom: NumpySnapshotExtension, toplogy):
    session = mn.io.MDAnalysisSession()
    remove_all_molecule_objects(session)
    universe = mda.Universe(
        data_dir / "martini" / toplogy, data_dir / "martini/pent/PENT2_100frames.xtc"
    )

    mol = session.show(universe, style="ribbon")

    pos_a = sample_attribute(mol, "position")
    bpy.context.scene.frame_set(3)
    pos_b = sample_attribute(mol, "position")

    assert not np.isclose(pos_a, pos_b).all()

    for att in mol.data.attributes.keys():
        assert snapshot_custom == sample_attribute(mol, att)

    for att in mol.data.attributes.keys():
        assert snapshot_custom == sample_attribute(mol, att)
