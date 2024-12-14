import bpy
import pytest
import MDAnalysis as mda
import numpy as np
from molecularnodes.entities.trajectory import dna
from .utils import NumpySnapshotExtension
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
            format=dna.OXDNAReader,
            topology_format=dna.OXDNAParser,
        )

    def test_read_as_universe(self, file_holl_top, file_holl_dat):
        u = mda.Universe(
            file_holl_top,
            file_holl_dat,
            format=dna.OXDNAReader,
            topology_format=dna.OXDNAParser,
        )
        assert u.atoms.n_atoms == 98

    def test_univ_as_traj(self, universe):
        traj = dna.OXDNA(universe)
        assert traj.universe
        assert not traj.object
        assert all([x in ["A", "C", "T", "G"] for x in traj.res_name])

    def test_univ_snapshot(self, universe: mda.Universe, snapshot_custom):
        traj = dna.OXDNA(universe)
        traj.create_object()
        for name in ["position", "res_name", "res_id", "chain_id"]:
            assert snapshot_custom == getattr(traj, name)

    def test_detect_new_top(self, file_top_old, file_top_new, file_top_new_custom):
        assert dna.OXDNAParser._is_new_topology(file_top_new)
        assert dna.OXDNAParser._is_new_topology(file_top_new_custom)
        assert not dna.OXDNAParser._is_new_topology(file_top_old)

    def test_topo_reading(
        self,
        file_top_old,
        file_top_new,
        file_top_new_custom,
    ):
        top_new = dna.OXDNAParser._read_topo_new(file_top_new)
        top_new_custom = dna.OXDNAParser._read_topo_new(file_top_new_custom)
        top_old = dna.OXDNAParser._read_topo_old(file_top_old)

        for top in [top_new, top_old, top_new_custom]:
            assert top.n_atoms == 12
            assert top.n_residues == 12

    @pytest.mark.parametrize("topfile", ["top_new", "top_new_custom", "top_old"])
    def test_comparing_topologies(self, snapshot, topfile, file_traj_old_new):
        u = mda.Universe(
            data_dir / f"oxdna/{topfile}.top",
            file_traj_old_new,
            topology_format=dna.OXDNAParser,
            format=dna.OXDNAReader,
        )
        traj = dna.OXDNA(u)
        traj.create_object()
        assert len(traj) == 12
        assert snapshot == traj.bonds.tolist()
        for att in ["res_id", "chain_id", "res_name"]:
            assert snapshot == traj.named_attribute(att).tolist()

    def test_reading_example(self):
        traj = dna.OXDNA(
            mda.Universe(
                data_dir / "CanDo2oxDNA/top.top",
                data_dir / "CanDo2oxDNA/traj.oxdna",
                topology_format=dna.OXDNAParser,
                format=dna.OXDNAReader,
            )
        )
        traj.create_object()
        assert len(np.unique(traj.named_attribute("res_id"))) == 15166
        assert len(np.unique(traj.named_attribute("chain_id"))) == 178

    def test_reload_lost_connection(self, snapshot, file_holl_top, file_holl_dat):
        u = mda.Universe(
            file_holl_top,
            file_holl_dat,
            topology_format=dna.OXDNAParser,
            format=dna.OXDNAReader,
        )
        traj = dna.OXDNA(u)
        traj.create_object()
        bpy.context.scene.frame_set(1)
        pos1 = traj.named_attribute("position")
        bpy.context.scene.frame_set(2)
        pos2 = traj.named_attribute("position")
        assert not np.allclose(pos1, pos2)

        del bpy.context.scene.MNSession.trajectories[traj.uuid]

        bpy.context.scene.frame_set(3)
        pos3 = traj.named_attribute("position")

        assert np.allclose(pos2, pos3)

        bpy.ops.mn.reload_trajectory()
        bpy.context.scene.frame_set(4)
        bpy.context.scene.frame_set(3)

        pos3 = traj.named_attribute("position")

        assert not np.allclose(pos2, pos3)
