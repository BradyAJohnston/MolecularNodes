import bpy
import os
import pytest
import MolecularNodes as mn
from MolecularNodes.mda import HAS_mda
if HAS_mda:
    import MDAnalysis as mda
import numpy as np
from .utils import get_verts, apply_mods, remove_all_molecule_objects

@pytest.mark.skipif(not HAS_mda, reason='MDAnalysis is not installed')
class MDATests:
    @pytest.fixture(scope='module')
    def mda_session():
        mda_session = mn.mda.MDAnalysisSession()
        return mda_session


    @pytest.fixture(scope='module')
    def universe():
        top = "tests/data/md_ppr/box.gro"
        traj = "tests/data/md_ppr/first_5_frames.xtc"
        u = mda.Universe(top, traj)
        return u


    def test_persistent_handlers_added(mda_session):
        assert bpy.app.handlers.load_post[-1].__name__ == '_rejuvenate_universe'
        assert bpy.app.handlers.save_pre[-1].__name__ == '_sync_universe'


    def test_create_mda_session(mda_session):
        assert mda_session is not None
        assert mda_session.uuid is not None
        assert mda_session.world_scale == 0.01


    def reload_mda_session(mda_session):
        with pytest.warns(UserWarning, match='The existing mda session'):
            mda_session_2 = mn.mda.create_session() 
        assert mda_session.uuid == mda_session_2.uuid


    def test_show_universe(snapshot, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe)
        obj = bpy.data.objects['atoms']
        verts = get_verts(obj, apply_modifiers = False)
        
        snapshot.assert_match(verts, 'md_gro_xtc_verts.txt')


    def test_same_name_atoms(snapshot, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe)

        with pytest.warns(UserWarning, match='The name of the object is changed'):
            mda_session.show(universe)
        obj_1 = bpy.data.objects['atoms']
        obj_2 = bpy.data.objects['atoms.001']
        verts_1 = get_verts(obj_1, apply_modifiers = False)
        verts_2 = get_verts(obj_2, apply_modifiers = False)
        
        assert(verts_1 == verts_2)


    def test_show_multiple_selection(snapshot, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        custom_selections = {'name_ca': 'name CA'}
        mda_session.show(universe,
                        name='protein',
                        selection='protein',
                        custom_selections=custom_selections)
        obj = bpy.data.objects['protein']
        verts = get_verts(obj, apply_modifiers = False)

        snapshot.assert_match(verts, 'md_gro_xtc_verts_protein.txt')
        obj_ca = bpy.data.objects['name_ca']
        verts_ca = get_verts(obj_ca, apply_modifiers = False)
        snapshot.assert_match(verts_ca, 'md_gro_xtc_verts_ca.txt')


    def test_trajectory_update(snapshot, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe)
        obj = bpy.data.objects['atoms']

        verts_frame_0 = get_verts(obj, apply_modifiers = False)
        snapshot.assert_match(verts_frame_0, 'md_gro_xtc_verts_frame_0.txt')

        # change blender frame to 1
        bpy.context.scene.frame_set(1)
        obj = bpy.data.objects['atoms']
        verts_frame_1 = get_verts(obj, apply_modifiers = False)
        snapshot.assert_match(verts_frame_1, 'md_gro_xtc_verts_frame_1.txt')

        assert(verts_frame_0 != verts_frame_1)


    def test_show_updated_atoms(snapshot, mda_session, universe):
        remove_all_molecule_objects(mda_session)
        updating_ag = universe.select_atoms('around 5 resid 1', updating=True)
        mda_session.show(updating_ag)

        obj = bpy.data.objects['atoms']

        verts_frame_0 = get_verts(obj, apply_modifiers = False)
        snapshot.assert_match(verts_frame_0, 'md_gro_xtc_verts_frame_0.txt')

        # change blender frame to 1
        bpy.context.scene.frame_set(1)
        print(mda_session.rep_names)
        obj = bpy.data.objects['atoms']
        verts_frame_1 = get_verts(obj, apply_modifiers = False)
        snapshot.assert_match(verts_frame_1, 'md_gro_xtc_verts_frame_1.txt')

        assert(verts_frame_0 != verts_frame_1)


    def test_save_persistance(snapshot, mda_session, universe, tmp_path):
        remove_all_molecule_objects(mda_session)
        mda_session.show(universe)
        # save
        bpy.ops.wm.save_as_mainfile(filepath=str(tmp_path / 'test.blend'))

        # reload
        remove_all_molecule_objects(mda_session)
        bpy.ops.wm.open_mainfile(filepath=str(tmp_path / 'test.blend'))
        obj = bpy.data.objects['atoms']
        verts_frame_0 = get_verts(obj, apply_modifiers = False)
        # change blender frame to 1
        bpy.context.scene.frame_set(1)
        obj = bpy.data.objects['atoms']
        verts_frame_1 = get_verts(obj, apply_modifiers = False)
        snapshot.assert_match(verts_frame_1, 'md_gro_xtc_verts_frame_1.txt')

        assert(verts_frame_0 != verts_frame_1)