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


@pytest.fixture()
def session():
    return mn.session.get_session()


@pytest.fixture()
def universe():
    topo = data_dir / "md_ppr/box.gro"
    traj = data_dir / "md_ppr/first_5_frames.xtc"
    return mda.Universe(topo, traj)


def test_add_trajectory(session, universe):
    # add Universe as trajectory
    t1 = session.add_trajectory(universe, name="u1")
    assert "u1" in bpy.data.objects
    assert t1._entity_type == mn.entities.base.EntityType.MD
    assert t1.object.mn.entity_type == t1._entity_type.value
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
    props = bpy.context.scene.mn
    assert len(props.entities) == 0
    t1 = session.add_trajectory(universe, name="u1")
    # test property addition
    assert len(props.entities) == 1
    entity = props.entities[0]
    # test property indexing
    assert entity.name == t1.uuid
    # verify entity type
    assert entity.type == "md"
    # test visibility changes
    # check initial visibility
    assert entity.visible
    assert bpy.data.objects["u1"].visible_get()
    entity.visible = False
    assert not bpy.data.objects["u1"].visible_get()
    # test property removal
    session.remove_trajectory("u1")
    assert len(props.entities) == 0
