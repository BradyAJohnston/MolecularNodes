import bpy
import os
import pytest
import molecularnodes as mn
from molecularnodes.blender.obj import get_attribute

import MDAnalysis as mda
import numpy as np
from .constants import data_dir
from .utils import sample_attribute, NumpySnapshotExtension

# mn.unregister()
# mn.register()


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

    @pytest.fixture(scope="module")
    def MNUniverse(self, universe):
        mnu = mn.io.md.MNUniverse(universe)
        mnu.create_model()
        return mnu

    @pytest.fixture(scope="module")
    def MNUniverse_with_bonds(self, universe_with_bonds):
        mnu = mn.io.md.MNUniverse(universe_with_bonds)
        mnu.create_model()
        return mnu

    @pytest.fixture(scope="module")
    def session(self):
        return bpy.context.scene.MNSession

    def test_include_bonds(self, MNUniverse_with_bonds):
        assert MNUniverse_with_bonds.object.data.edges.items() != []

    def test_attributes_added(self, MNUniverse):
        attributes = MNUniverse.object.data.attributes.keys()
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

    def test_trajectory_update(self, snapshot_custom, MNUniverse):
        bob = MNUniverse.object
        print(f"{bpy.context.scene.MNSession.universes=}")

        bpy.context.scene.frame_set(0)
        pos_a = get_attribute(bob, "position")
        assert snapshot_custom == pos_a

        bpy.context.scene.frame_set(4)
        bob.data.update()
        pos_b = get_attribute(bob, "position")
        assert snapshot_custom == pos_b

        assert not np.allclose(pos_a, pos_b)

    @pytest.mark.parametrize("interpolate", [True, False])
    def test_subframes(self, MNUniverse, interpolate):
        bob = MNUniverse.object
        bpy.context.scene.frame_set(0)
        verts_a = get_attribute(bob, "position")

        bpy.context.scene.frame_set(1)
        verts_b = get_attribute(bob, "position")
        # should be no difference because not using subframes
        assert not np.isclose(verts_a, verts_b).all()

        for subframes in [1, 2, 3, 4]:
            frame = 1
            fraction = frame % (subframes + 1) / (subframes + 1)
            bob.mn.subframes = subframes
            bob.mn.interpolate = interpolate

            bpy.context.scene.frame_set(frame)
            verts_c = get_attribute(bob, "position")

            if interpolate:
                # now using subframes, there should be a difference
                assert not np.allclose(verts_b, verts_c)
                assert np.allclose(verts_c, mn.utils.lerp(verts_a, verts_b, t=fraction))
            else:
                assert np.allclose(verts_b, verts_c)

    def test_save_persistance(
        self,
        snapshot_custom: NumpySnapshotExtension,
        tmp_path,
        universe,
        session: mn.session.MNSession,
    ):
        session.clear()
        mnu = mn.io.MNUniverse(universe)
        mnu.create_model()
        object_name = mnu.object.name
        bpy.context.scene.frame_set(0)
        filepath = str(tmp_path / "test.blend")

        # test that we can save the file and it is created only after saving
        assert not os.path.exists(session.stashpath(filepath))
        bpy.ops.wm.save_as_mainfile(filepath=filepath)
        assert os.path.exists(filepath)
        assert os.path.exists(session.stashpath(filepath))
        bpy.ops.wm.open_mainfile(filepath=filepath)

        bob = bpy.data.objects[object_name]
        verts_frame_0 = mn.blender.obj.get_attribute(bob, "position")
        bpy.context.scene.frame_set(4)
        verts_frame_4 = mn.blender.obj.get_attribute(bob, "position")

        assert snapshot_custom == verts_frame_4
        assert not np.allclose(verts_frame_0, verts_frame_4)


@pytest.mark.parametrize("toplogy", ["pent/prot_ion.tpr", "pent/TOPOL2.pdb"])
def test_martini(snapshot_custom: NumpySnapshotExtension, toplogy):
    universe = mda.Universe(
        data_dir / "martini" / toplogy, data_dir / "martini/pent/PENT2_100frames.xtc"
    )
    mnu = mn.io.MNUniverse(universe)
    mnu.create_model()
    bob = mnu.object
    bpy.context.scene.frame_set(0)
    pos_a = get_attribute(bob, "position")

    bpy.context.scene.frame_set(50)
    pos_b = get_attribute(bob, "position")
    assert not np.allclose(pos_a, pos_b)

    for att in bob.data.attributes.keys():
        assert snapshot_custom == sample_attribute(bob, att)

    for att in bob.data.attributes.keys():
        assert snapshot_custom == sample_attribute(bob, att)
