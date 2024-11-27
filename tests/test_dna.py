import numpy as np
import molecularnodes as mn
import pytest
import MDAnalysis as mda
from molecularnodes.entities.trajectory import dna
from .utils import NumpySnapshotExtension
from .constants import data_dir


class TestOXDNAReading:
    @pytest.fixture(scope="module")
    def filepath_top(self):
        return data_dir / "oxdna/holliday.top"

    @pytest.fixture(scope="module")
    def filepath_traj_dat(self):
        return data_dir / "oxdna/holliday.dat"

    @pytest.fixture(scope="module")
    def universe(self, filepath_top, filepath_traj_dat):
        return mda.Universe(
            filepath_top,
            filepath_traj_dat,
            format=dna.OXDNAReader,
            topology_format=dna.OXDNAParser,
        )

    def test_read_as_universe(self, filepath_top, filepath_traj_dat):
        u = mda.Universe(
            filepath_top,
            filepath_traj_dat,
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
