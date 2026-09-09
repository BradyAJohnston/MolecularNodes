import bpy
import MDAnalysis as mda
import numpy as np
import pytest
from databpy.object import LinkedObjectError
import molecularnodes as mn
from molecularnodes.entities.molecule import oxdna
from .constants import data_dir

pytestmark = pytest.mark.filterwarnings(
    "ignore:.*no reference attributes.*:UserWarning"
)


class TestOXDNAReading:
    # def filepath(self, file):
    #     return data_dir / f"oxdna/{file}"

    # @pytest.fixture(
    #     scope="module",
    #     params=[
    #         "linear.top",
    #         "linear_custom.top",
    #         "linear_traj.dat",
    #         "linear_old.top",
    #         "linear_old_traj.dat",
    #         "holliday_old.top",
    #         "holliday_old_traj.dat",
    #         "minicircle.top",
    #         "minicircle.dat",
    #         "minicircle_old.top",
    #         "minicircle_old.dat",
    #         "origami_old.top",
    #         "origami_old.dat",
    #     ],
    # )
    # def file_name(self, request):
    #     return request.param

    # def test_read_as_universe(self, snapshot, file_holl_top_old, file_holl_traj_old):
    #     u = mda.Universe(
    #         file_holl_top_old,
    #         file_holl_traj_old,
    #         format=oxdna.OXDNAReader,
    #         topology_format=oxdna.OXDNAParser,
    #     )
    #     assert snapshot == u.atoms.n_atoms

    @pytest.fixture(scope="module")
    def file_lin_top(self):
        return data_dir / "oxdna/linear.top"

    @pytest.fixture(scope="module")
    def file_lin_top_custom(self):
        return data_dir / "oxdna/linear_custom.top"

    @pytest.fixture(scope="module")
    def file_lin_traj(self):
        return data_dir / "oxdna/linear_traj.dat"

    @pytest.fixture(scope="module")
    def file_lin_top_old(self):
        return data_dir / "oxdna/linear_old.top"

    @pytest.fixture(scope="module")
    def file_lin_traj_old(self):
        return data_dir / "oxdna/linear_old_traj.dat"

    @pytest.fixture(scope="module")
    def file_holl_top_old(self):
        return data_dir / "oxdna/holliday_old.top"

    @pytest.fixture(scope="module")
    def file_holl_traj_old(self):
        return data_dir / "oxdna/holliday_old_traj.dat"

    @pytest.fixture(scope="module")
    def file_circ_top(self):
        return data_dir / "oxdna/minicircle.top"

    @pytest.fixture(scope="module")
    def file_circ_conf(self):
        return data_dir / "oxdna/minicircle.dat"

    @pytest.fixture(scope="module")
    def file_circ_top_old(self):
        return data_dir / "oxdna/minicircle_old.top"

    @pytest.fixture(scope="module")
    def file_circ_conf_old(self):
        return data_dir / "oxdna/minicircle_old.dat"

    @pytest.fixture(scope="module")
    def file_origami_top_old(self):
        return data_dir / "oxdna/origami_old.top"

    @pytest.fixture(scope="module")
    def file_origami_conf_old(self):
        return data_dir / "oxdna/origami_old.dat"

    def test_read_as_universe(self, snapshot, file_holl_top_old, file_holl_traj_old):
        u = mda.Universe(
            file_holl_top_old,
            file_holl_traj_old,
            format=oxdna.OXDNAReader,
            topology_format=oxdna.OXDNAParser,
        )
        assert snapshot == u.atoms.n_atoms

    def test_univ_as_traj(self, universe):
        traj = oxdna.OXDNA(universe, create_object=False)
        assert traj.universe
        with pytest.raises(LinkedObjectError):
            traj.object
        assert all([x in ["A", "C", "T", "G"] for x in traj.atoms.resnames])

    def test_univ_snapshot(self, universe: mda.Universe, snapshot_custom):
        traj = oxdna.OXDNA(universe)
        for name in ["position", "res_name", "res_id", "chain_id"]:
            assert snapshot_custom == str(traj[name])

    def test_detect_new_top(self, file_lin_top_old, file_lin_top, file_lin_top_custom):
        assert oxdna.OXDNAParser._is_new_topology(file_lin_top)
        assert oxdna.OXDNAParser._is_new_topology(file_lin_top_custom)
        assert not oxdna.OXDNAParser._is_new_topology(file_lin_top_old)

    def test_topo_reading(
        self, snapshot, file_lin_top_old, file_lin_top, file_lin_top_custom
    ):
        top_new = oxdna.OXDNAParser._read_topo_new(file_lin_top)
        top_new_custom = oxdna.OXDNAParser._read_topo_new(file_lin_top_custom)
        top_old = oxdna.OXDNAParser._read_topo_old(file_lin_top_old)

        for top in [top_new, top_old, top_new_custom]:
            assert snapshot == top.n_atoms
            assert snapshot == top.n_residues

    @pytest.mark.parametrize(
        "topfile, trajfile",
        [
            ("linear", "linear_traj"),
            ("linear_custom", "linear_traj"),
            ("linear_old", "linear_old_traj"),
        ],
    )
    def test_comparing_topologies(self, snapshot, topfile, trajfile):
        u = mda.Universe(
            data_dir / f"oxdna/{topfile}.top",
            data_dir / f"oxdna/{trajfile}.dat",
            topology_format=oxdna.OXDNAParser,
            format=oxdna.OXDNAReader,
        )
        traj = oxdna.OXDNA(u)
        assert snapshot == len(traj)
        assert snapshot == traj.atoms.bonds.indices.tolist()
        for att in ["res_id", "chain_id", "res_name"]:
            assert snapshot == str(traj[att])

    @pytest.mark.parametrize(
        "topfile, trajfile",
        [
            ("minicircle", "minicircle"),
            ("minicircle_old", "minicircle_old"),
            ("linear", "linear_traj"),
            ("linear_old", "linear_old_traj"),
        ],
    )
    def test_reading_ligation(self, snapshot, topfile, trajfile):
        u = mda.Universe(
            data_dir / f"oxdna/{topfile}.top",
            data_dir / f"oxdna/{trajfile}.dat",
            topology_format=oxdna.OXDNAParser,
            format=oxdna.OXDNAReader,
        )
        traj = oxdna.OXDNA(u)
        assert snapshot == traj.atoms.bonds.indices.tolist()

    @pytest.mark.parametrize("topfile, trajfile", [("origami_old", "origami_old")])
    def test_reading_example(self, snapshot, topfile, trajfile):
        u = mda.Universe(
            data_dir / f"oxdna/{topfile}.top",
            data_dir / f"oxdna/{trajfile}.dat",
            topology_format=oxdna.OXDNAParser,
            format=oxdna.OXDNAReader,
        )
        traj = oxdna.OXDNA(u)
        for att in ["res_id", "chain_id"]:
            assert snapshot == len(np.unique(traj.named_attribute[att]))

    def test_session_register(self, file_holl_top_old, file_holl_traj_old):
        session = mn.session.get_session()
        u = mda.Universe(
            file_holl_top_old,
            file_holl_traj_old,
            topology_format=oxdna.OXDNAParser,
            format=oxdna.OXDNAReader,
        )
        traj = oxdna.OXDNA(u)

        assert isinstance(session.get(traj.uuid), oxdna.OXDNA)
        assert traj._mn_entity_type == mn.entities.base.EntityType.MD_OXDNA.value

    def test_reload_lost_connection(self, file_holl_top_old, file_holl_traj_old):
        session = mn.session.get_session()
        u = mda.Universe(
            file_holl_top_old,
            file_holl_traj_old,
            topology_format=oxdna.OXDNAParser,
            format=oxdna.OXDNAReader,
        )
        traj = oxdna.OXDNA(u)
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
