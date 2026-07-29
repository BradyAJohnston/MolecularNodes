import bpy
import MDAnalysis as mda
import pytest
import molecularnodes as mn
from .constants import data_dir


def test_session_present():
    assert isinstance(mn.session.get_session(), mn.session.MNSession)


def test_persistent_handlers_added():
    load_handlers = [handler.__name__ for handler in bpy.app.handlers.load_post]
    save_handlers = [handler.__name__ for handler in bpy.app.handlers.save_post]
    assert "_pickle" in save_handlers
    assert "_load" in load_handlers


def test_entity_registered():
    session = mn.session.get_session()
    assert len(session.entities) == 0

    mol = mn.Molecule.fetch("1BNA", cache=data_dir)

    assert mol.uuid in session.entities
    assert isinstance(session.get(mol.uuid), mn.Molecule)
    assert len(session.entities) == 1


def test_reload_molecule_from_file(tmp_path):
    from molecularnodes.entities.reload import can_reload, reload_entity

    session = mn.session.get_session()
    mol = mn.Molecule.load(data_dir / "1cd3.cif")
    obj = mol.object
    # the source file is recorded so the entity can be reloaded
    assert obj.mn.filepath != ""
    assert can_reload(obj)

    # simulate a lost session link (e.g. re-opening a .blend without the pickle)
    session.remove_entity(obj.uuid)
    assert session.get(obj.uuid) is None

    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.Molecule)
    assert session.get(obj.uuid) is reloaded
    assert reloaded.object is obj


def test_reload_molecule_from_code():
    from molecularnodes.entities.reload import reload_entity

    session = mn.session.get_session()
    mol = mn.Molecule.fetch("1BNA", cache=data_dir)
    obj = mol.object
    assert obj.mn.code == "1BNA"

    session.remove_entity(obj.uuid)
    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.Molecule)
    assert session.get(obj.uuid) is reloaded


def test_reload_trajectory_from_file():
    from molecularnodes.entities.reload import can_reload, reload_entity

    session = mn.session.get_session()
    traj = mn.entities.trajectory.load(
        top=data_dir / "md_ppr/box.gro",
        traj=data_dir / "md_ppr/first_5_frames.xtc",
    )
    obj = traj.object
    assert obj.mn.filepath_topology != ""
    assert can_reload(obj)

    session.remove_entity(obj.uuid)
    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.Trajectory)
    assert session.get(obj.uuid) is reloaded


def test_reload_density_from_file():
    from molecularnodes.entities.reload import can_reload, reload_entity

    session = mn.session.get_session()
    density = mn.entities.density.load(data_dir / "emd_24805.map.gz")
    obj = density.object
    assert obj.mn.filepath != ""
    assert can_reload(obj)

    session.remove_entity(obj.uuid)
    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.entities.density.Density)
    assert session.get(obj.uuid) is reloaded


def test_reload_ensemble_star_from_file():
    from molecularnodes.entities.reload import can_reload, reload_entity

    session = mn.session.get_session()
    ensemble = mn.entities.ensemble.StarFile.load(data_dir / "starfile/relion.star")
    obj = ensemble.object
    assert can_reload(obj)

    session.remove_entity(obj.uuid)
    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.entities.ensemble.StarFile)
    assert session.get(obj.uuid) is reloaded


def test_reload_ensemble_cellpack_from_file():
    from molecularnodes.entities.reload import can_reload, reload_entity

    session = mn.session.get_session()
    ensemble = mn.entities.ensemble.CellPack.load(
        file_path=data_dir / "cellpack/square1.bcif",
        name="CellPack",
        node_setup=False,
    )
    obj = ensemble.object
    assert obj.mn.filepath != ""
    assert can_reload(obj)

    session.remove_entity(obj.uuid)
    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.entities.ensemble.CellPack)
    assert session.get(obj.uuid) is reloaded


@pytest.fixture()
def session():
    return mn.session.get_session()


@pytest.fixture()
def universe():
    topo = data_dir / "md_ppr/box.gro"
    traj = data_dir / "md_ppr/first_5_frames.xtc"
    return mda.Universe(topo, traj)


@pytest.mark.filterwarnings("ignore:.*Empty string to select atoms.*:UserWarning")
def test_add_trajectory(session, universe):
    # add Universe as trajectory
    t1 = session.add_trajectory(universe, name="u1")
    assert "u1" in bpy.data.objects
    assert t1._mn_entity_type == mn.entities.base.EntityType.MD.value
    # add AtomGroup as trajectory
    ag = universe.select_atoms("name CA")
    session.add_trajectory(ag, name="ag1")
    assert "ag1" in bpy.data.objects


def test_remove_trajectory(session, universe):
    t1 = session.add_trajectory(universe, name="u1")
    assert "u1" in bpy.data.objects
    # remove by trajectory instance
    session.remove_trajectory(t1)
    assert "u1" not in bpy.data.objects
    session.add_trajectory(universe, name="u2")
    assert "u2" in bpy.data.objects
    # remove by trajectory name
    session.remove_trajectory("u2")
    t3 = session.add_trajectory(universe, name="u3")
    assert "u3" in bpy.data.objects
    # remove trajectory from UI using operator
    bpy.ops.mn.session_remove_item("EXEC_DEFAULT", uuid=t3.uuid)
    assert "u3" not in bpy.data.objects


def test_get_trajectory(session, universe):
    t1 = session.add_trajectory(universe, name="u1")
    t2 = session.get_trajectory("u1")
    assert t1 == t2


def test_entity_blender_properties(session, universe):
    assert len(session.entities) == 0
    t1 = session.add_trajectory(universe, name="u1")
    # entity is tracked in the session, keyed by uuid
    assert len(session.entities) == 1
    assert t1.uuid in session.entities
    obj = bpy.data.objects["u1"]
    # entity type is stored on the object properties
    assert obj.mn.entity_type == "md"
    assert obj.uuid == t1.uuid
    # visibility is driven through the object properties
    assert obj.mn.visible
    assert obj.visible_get()
    obj.mn.visible = False
    assert not obj.visible_get()
    # removal drops the entity from the session
    session.remove_trajectory("u1")
    assert len(session.entities) == 0
