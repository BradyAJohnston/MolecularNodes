import bpy
import MDAnalysis as mda
import numpy as np
import pytest
from databpy.object import LinkedObjectError
import molecularnodes as mn
from molecularnodes.entities.trajectory import oxdna
from .constants import data_dir


class TestOXDNAReading:
    @pytest.fixture(scope="module")
    def file_holl_top(self):
        return data_dir / "oxdna/holliday.top"

    @pytest.fixture(scope="module")
    def file_top_new(self):
        return data_dir / "oxdna/top_new.top"

    @pytest.fixture(scope="module")
    def file_top_new_custom(self):
        return data_dir / "oxdna/top_new_custom.top"

    @pytest.fixture(scope="module")
    def file_top_old(self):
        return data_dir / "oxdna/top_old.top"

    @pytest.fixture(scope="module")
    def file_traj_old_new(self):
        return data_dir / "oxdna/traj_old_new.dat"

    @pytest.fixture(scope="module")
    def file_holl_dat(self):
        return data_dir / "oxdna/holliday.dat"

    @pytest.fixture(scope="module")
    def universe(self, file_holl_top, file_holl_dat):
        return mda.Universe(
            file_holl_top,
            file_holl_dat,
            format=oxdna.OXDNAReader,
            topology_format=oxdna.OXDNAParser,
        )

    def test_read_as_universe(self, file_holl_top, file_holl_dat):
        u = mda.Universe(
            file_holl_top,
            file_holl_dat,
            format=oxdna.OXDNAReader,
            topology_format=oxdna.OXDNAParser,
        )
        assert u.atoms.n_atoms == 98

    def test_univ_as_traj(self, universe):
        traj = oxdna.OXDNA(universe)
        assert traj.universe
        with pytest.raises(LinkedObjectError):
            traj.object
        assert all([x in ["A", "C", "T", "G"] for x in traj.res_name])

    def test_univ_snapshot(self, universe: mda.Universe, snapshot_custom):
        traj = oxdna.OXDNA(universe)
        traj.create_object()
        for name in ["position", "res_name", "res_id", "chain_id"]:
            assert snapshot_custom == getattr(traj, name)

    def test_detect_new_top(self, file_top_old, file_top_new, file_top_new_custom):
        assert oxdna.OXDNAParser._is_new_topology(file_top_new)
        assert oxdna.OXDNAParser._is_new_topology(file_top_new_custom)
        assert not oxdna.OXDNAParser._is_new_topology(file_top_old)

    def test_topo_reading(
        self,
        file_top_old,
        file_top_new,
        file_top_new_custom,
    ):
        top_new = oxdna.OXDNAParser._read_topo_new(file_top_new)
        top_new_custom = oxdna.OXDNAParser._read_topo_new(file_top_new_custom)
        top_old = oxdna.OXDNAParser._read_topo_old(file_top_old)

        for top in [top_new, top_old, top_new_custom]:
            assert top.n_atoms == 12
            assert top.n_residues == 12

    @pytest.mark.parametrize("topfile", ["top_new", "top_new_custom", "top_old"])
    def test_comparing_topologies(self, snapshot, topfile, file_traj_old_new):
        u = mda.Universe(
            data_dir / f"oxdna/{topfile}.top",
            file_traj_old_new,
            topology_format=oxdna.OXDNAParser,
            format=oxdna.OXDNAReader,
        )
        traj = oxdna.OXDNA(u)
        traj.create_object()
        assert len(traj) == 12
        assert snapshot == traj.bonds.tolist()
        for att in ["res_id", "chain_id", "res_name"]:
            assert snapshot == traj.named_attribute(att).tolist()

    def test_reading_example(self):
        traj = oxdna.OXDNA(
            mda.Universe(
                data_dir / "CanDo2oxDNA/top.top",
                data_dir / "CanDo2oxDNA/traj.oxdna",
                topology_format=oxdna.OXDNAParser,
                format=oxdna.OXDNAReader,
            )
        )
        traj.create_object()
        assert len(np.unique(traj.named_attribute("res_id"))) == 15166
        assert len(np.unique(traj.named_attribute("chain_id"))) == 178

    def test_session_register(self, file_holl_top, file_holl_dat):
        session = mn.session.get_session()
        u = mda.Universe(
            file_holl_top,
            file_holl_dat,
            topology_format=oxdna.OXDNAParser,
            format=oxdna.OXDNAReader,
        )
        traj = oxdna.OXDNA(u)
        traj.create_object()

        assert isinstance(session.get(traj.uuid), oxdna.OXDNA)

    def test_reload_lost_connection(self, snapshot, file_holl_top, file_holl_dat):
        session = mn.session.get_session()
        u = mda.Universe(
            file_holl_top,
            file_holl_dat,
            topology_format=oxdna.OXDNAParser,
            format=oxdna.OXDNAReader,
        )
        traj = oxdna.OXDNA(u)
        traj.create_object()
        obj_name = traj.name
        bpy.context.scene.frame_set(1)
        pos1 = traj.position
        bpy.context.scene.frame_set(2)
        pos2 = traj.position
        assert not np.allclose(pos1, pos2)

        traj_old = session.entities.pop(traj.uuid)

        # the position shouldn't change as we have removed the traj from the session
        bpy.context.scene.frame_set(3)
        pos3 = traj.position
        assert np.allclose(pos2, pos3)
        del traj

        bpy.data.objects[obj_name].select_set(True)
        bpy.ops.mn.reload_trajectory()

        # when reloading the object, a brand new traj had to be created, which updates
        # the uuid on the object, so the old traj will not longer be able to find any
        # matching object and instead we'll have to look back up a new traj based on the
        # the object's uuid
        with pytest.raises(LinkedObjectError):
            traj_old.object.name

        traj = bpy.context.scene.MNSession.get(bpy.context.active_object.uuid)
        assert traj is not None

        pos3 = traj.position

        assert not np.allclose(pos2, pos3)
