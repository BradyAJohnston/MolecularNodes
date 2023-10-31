import bpy
import os
import pytest
import molecularnodes as mn
from molecularnodes.mda import HAS_mda

if HAS_mda:
    import MDAnalysis as mda
import numpy as np
from .constants import (
    test_data_directory
)
from .utils import get_verts, apply_mods, remove_all_molecule_objects

@pytest.mark.skipif(not HAS_mda, reason="MDAnalysis is not installed")
class TestMDA:
    @pytest.fixture(scope="module")
    def mda_session(self):
        mda_session = mn.mda.MDAnalysisSession()
        return mda_session

    @pytest.fixture(scope="module")
    def universe(self):
        top = test_data_directory / "md_ppr/box.gro"
        traj = test_data_directory / "md_ppr/first_5_frames.xtc"
        u = mda.Universe(top, traj)
        return u

    @pytest.fixture(scope="module")
    def universe_with_bonds(self):
        top = test_data_directory / "md_ppr/md.tpr"
        traj = test_data_directory / "md_ppr/md.gro"
        u = mda.Universe(top, traj)
        return u

    def test_persistent_handlers_added(self, mda_session):
        assert bpy.app.handlers.load_post[-1].__name__ == "_rejuvenate_universe"
        assert bpy.app.handlers.save_pre[-1].__name__ == "_sync_universe"

    def test_create_mda_session(self, mda_session):
        assert mda_session is not None
        assert mda_session.uuid is not None
        assert mda_session.world_scale == 0.01

    def reload_mda_session(self, mda_session):
        with pytest.warns(UserWarning, match="The existing mda session"):
            mda_session_2 = mn.mda.create_session()
        assert mda_session.uuid == mda_session_2.uuid

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_show_universe(self, snapshot, in_memory, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe, in_memory=in_memory)
        obj = bpy.data.objects["atoms"]
        verts = get_verts(obj, apply_modifiers=False)

        snapshot.assert_match(verts, "md_gro_xtc_verts.txt")

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_same_name_atoms(self, snapshot, in_memory, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe, in_memory=in_memory)

        with pytest.warns(UserWarning, match="The name of the object is changed"):
            mda_session.show(universe, in_memory=in_memory)

        obj_1 = bpy.data.objects["atoms"]
        obj_2 = bpy.data.objects["atoms.001"]
        verts_1 = get_verts(obj_1, apply_modifiers=False)
        verts_2 = get_verts(obj_2, apply_modifiers=False)

        assert verts_1 == verts_2

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_show_multiple_selection(self, snapshot, in_memory, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        custom_selections = {"name_ca": "name CA"}
        mda_session.show(
            universe,
            in_memory=in_memory,
            name="protein",
            selection="protein",
            custom_selections=custom_selections,
        )
        obj = bpy.data.objects["protein"]
        verts = get_verts(obj, apply_modifiers=False)

        snapshot.assert_match(verts, "md_gro_xtc_verts_protein.txt")

        # different bahavior in_memory or not.
        if not in_memory:
            obj_ca = bpy.data.objects["name_ca"]
            verts_ca = get_verts(obj_ca, apply_modifiers=False)
            snapshot.assert_match(verts_ca, "md_gro_xtc_verts_ca.txt")
        else:
            # attribute is added as name_ca.
            assert "name_ca" in obj.data.attributes.keys()

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_include_bonds(self, in_memory, mda_session, universe_with_bonds):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe_with_bonds, in_memory=in_memory, include_bonds=False)
        obj = bpy.data.objects["atoms"]
        assert obj.data.edges.items() == []

        remove_all_molecule_objects(mda_session)
        mda_session.show(universe_with_bonds, in_memory=in_memory, include_bonds=True)
        obj = bpy.data.objects["atoms"]
        assert obj.data.edges.items() != []

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_attributes_added(self, in_memory, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe, in_memory=in_memory, include_bonds=False)
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
    def test_trajectory_update(self, snapshot, in_memory, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe, in_memory=in_memory)
        obj = bpy.data.objects["atoms"]

        verts_frame_0 = get_verts(obj, apply_modifiers=False)
        snapshot.assert_match(verts_frame_0, "md_gro_xtc_verts_frame_0.txt")

        # change blender frame to 1
        bpy.context.scene.frame_set(1)
        obj = bpy.data.objects["atoms"]
        verts_frame_1 = get_verts(obj, apply_modifiers=False)
        snapshot.assert_match(verts_frame_1, "md_gro_xtc_verts_frame_1.txt")

        assert verts_frame_0 != verts_frame_1

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_show_updated_atoms(self, snapshot, in_memory, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        updating_ag = universe.select_atoms("around 5 resid 1", updating=True)
        mda_session.show(updating_ag, in_memory=in_memory)

        obj = bpy.data.objects["atoms"]

        verts_frame_0 = get_verts(obj, apply_modifiers=False)
        snapshot.assert_match(verts_frame_0, "md_gro_xtc_verts_frame_0.txt")

        # change blender frame to 1
        bpy.context.scene.frame_set(1)
        print(mda_session.rep_names)
        obj = bpy.data.objects["atoms"]
        verts_frame_1 = get_verts(obj, apply_modifiers=False)
        snapshot.assert_match(verts_frame_1, "md_gro_xtc_verts_frame_1.txt")

        assert verts_frame_0 != verts_frame_1

    @pytest.mark.parametrize("in_memory", [False, True])
    def test_update_deleted_objects(self, snapshot, in_memory, mda_session, universe):
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
        self, snapshot, tmp_path, in_memory, mda_session, universe
    ):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe, in_memory=in_memory)
        # save
        bpy.ops.wm.save_as_mainfile(filepath=str(tmp_path / "test.blend"))

        assert os.path.exists(f"{mda_session.session_tmp_dir}/{mda_session.uuid}.pkl")

        # reload
        remove_all_molecule_objects(mda_session)
        bpy.ops.wm.open_mainfile(filepath=str(tmp_path / "test.blend"))
        obj = bpy.data.objects["atoms"]
        verts_frame_0 = get_verts(obj, apply_modifiers=False)
        # change blender frame to 1
        bpy.context.scene.frame_set(1)
        obj = bpy.data.objects["atoms"]
        verts_frame_1 = get_verts(obj, apply_modifiers=False)
        snapshot.assert_match(verts_frame_1, "md_gro_xtc_verts_frame_1.txt")

        assert verts_frame_0 != verts_frame_1
