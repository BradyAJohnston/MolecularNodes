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
    def filepath_top_new(self):
        return data_dir / "oxdna/top_new.top"

    @pytest.fixture(scope="module")
    def filepath_top_new_custom(self):
        return data_dir / "oxdna/top_new_custom.top"

    @pytest.fixture(scope="module")
    def filepath_top_old(self):
        return data_dir / "oxdna/top_old.top"

    @pytest.fixture(scope="module")
    def filepath_traj_old_new(self):
        return data_dir / "oxdna/traj_old_new.dat"

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

    def test_detect_new_top(
        self, filepath_top_old, filepath_top_new, filepath_top_new_custom
    ):
        assert dna.OXDNAParser._is_new_topology(filepath_top_new)
        assert dna.OXDNAParser._is_new_topology(filepath_top_new_custom)
        assert not dna.OXDNAParser._is_new_topology(filepath_top_old)

    def test_topo_reading(
        self,
        filepath_top_old,
        filepath_top_new,
        filepath_top_new_custom,
    ):
        top_new = dna.OXDNAParser._read_topo_new(filepath_top_new)
        top_new_custom = dna.OXDNAParser._read_topo_new(filepath_top_new_custom)
        top_old = dna.OXDNAParser._read_topo_old(filepath_top_old)

        for top in [top_new, top_old, top_new_custom]:
            assert top.n_atoms == 12
            assert top.n_residues == 12

    @pytest.mark.parametrize("topfile", ["top_new", "top_new_custom", "top_old"])
    def test_comparing_topologies(self, snapshot, topfile, filepath_traj_old_new):
        u = mda.Universe(
            data_dir / f"oxdna/{topfile}.top",
            filepath_traj_old_new,
            topology_format=dna.OXDNAParser,
            format=dna.OXDNAReader,
        )
        traj = dna.OXDNA(u)
        traj.create_object()
        assert len(traj) == 12
        assert snapshot == traj.bonds.tolist()
        for att in ["res_id", "chain_id", "res_name"]:
            assert snapshot == traj.named_attribute(att).tolist()
