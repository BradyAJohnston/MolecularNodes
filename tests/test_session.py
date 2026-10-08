import bpy
import MDAnalysis as mda
import pytest
import molecularnodes as mn
from molecularnodes.session import MNSession
from .constants import data_dir


def test_session_present():
    assert isinstance(mn.session.get_session(), mn.session.MNSession)


def test_persistent_handlers_added():
    load_handlers = [handler.__name__ for handler in bpy.app.handlers.load_post]
    save_handlers = [handler.__name__ for handler in bpy.app.handlers.save_post]
    assert "_pickle" in save_handlers
    assert "_load" in load_handlers


def test_entity_registered(fetch):
    session = mn.session.get_session()
    assert len(session.entities) == 0

    mol = fetch("1BNA")

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


def test_reload_molecule_from_code(fetch, monkeypatch):
    import functools
    from molecularnodes.download import StructureDownloader
    from molecularnodes.entities import reload
    from molecularnodes.entities.reload import reload_entity

    # reload re-fetches the code into the default cache; read tests/data instead
    monkeypatch.setattr(
        reload,
        "StructureDownloader",
        functools.partial(StructureDownloader, cache=data_dir),
    )

    session = mn.session.get_session()
    mol = fetch("1BNA")
    obj = mol.object
    assert obj.mn.code == "1BNA"

    session.remove_entity(obj.uuid)
    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.Molecule)
    assert session.get(obj.uuid) is reloaded


def test_reload_trajectory_from_file():
    from molecularnodes.entities.reload import can_reload, reload_entity

    session = mn.session.get_session()
    traj = mn.Molecule.load(
        topology=data_dir / "md_ppr/box.gro",
        coordinates=data_dir / "md_ppr/first_5_frames.xtc",
    )
    obj = traj.object
    assert obj.mn.filepath_topology != ""
    assert can_reload(obj)

    session.remove_entity(obj.uuid)
    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.Molecule)
    assert session.get(obj.uuid) is reloaded


def test_reload_density_from_file():
    from molecularnodes.entities.reload import can_reload, reload_entity

    session = mn.session.get_session()
    density = mn.entities.density.Grids.load(data_dir / "emd_24805.map.gz")
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


def test_session_pickle_roundtrip_cellpack(tmp_path):
    session = mn.session.get_session()
    ensemble = mn.entities.ensemble.CellPack.load(
        data_dir / "cellpack/square1.bcif", node_setup=False
    )
    blend_path = tmp_path / "test.blend"

    session.pickle(blend_path)

    session.clear()
    session.load(blend_path)
    restored = session.get(ensemble.uuid)
    assert isinstance(restored, mn.entities.ensemble.CellPack)
    assert restored.name == ensemble.name
    assert restored.instance_collection is not None


def test_session_pickle_skips_unpicklable_entity(tmp_path):
    import threading

    session = mn.session.get_session()
    good = mn.Molecule.load(data_dir / "1cd3.cif")
    bad = mn.Molecule.load(data_dir / "1cd3.cif")
    bad._unpicklable = threading.Lock()
    blend_path = tmp_path / "test.blend"

    with pytest.warns(UserWarning, match="cannot be serialized"):
        session.pickle(blend_path)

    # the failed entity is skipped, but doesn't stop the rest of the session saving
    session.clear()
    session.load(blend_path)
    assert session.get(good.uuid) is not None
    assert session.get(bad.uuid) is None


@pytest.fixture()
def session():
    return mn.session.get_session()


@pytest.fixture()
def universe():
    topo = data_dir / "md_ppr/box.gro"
    traj = data_dir / "md_ppr/first_5_frames.xtc"
    return mda.Universe(topo, traj)


def test_multi_file_trajectory_paths_relativized():
    from contextlib import chdir
    from pathlib import Path
    from molecularnodes.session import (
        _make_trajectory_paths_absolute,
        _make_trajectory_paths_relative,
    )

    topo = data_dir / "md_ppr/box.gro"
    coords = data_dir / "md_ppr/first_5_frames.xtc"
    traj = mn.Molecule(mda.Universe(topo, [coords, coords]), name="multi")
    n_frames = traj.universe.trajectory.n_frames
    trajectories = {traj.uuid: traj}

    with chdir(data_dir):
        _make_trajectory_paths_relative(trajectories)
        filenames = [Path(f) for f in traj.universe.trajectory.filenames]
        assert len(filenames) == 2
        assert all(not f.is_absolute() for f in filenames)
        assert traj.universe.trajectory.n_frames == n_frames

        _make_trajectory_paths_absolute(trajectories)
        filenames = [Path(f) for f in traj.universe.trajectory.filenames]
        assert len(filenames) == 2
        assert all(f.is_absolute() for f in filenames)
        assert traj.universe.trajectory.n_frames == n_frames


def test_session_pickle_roundtrip_multi_file_trajectory(tmp_path):
    from contextlib import chdir
    from pathlib import Path

    topo = data_dir / "md_ppr/box.gro"
    coords = data_dir / "md_ppr/first_5_frames.xtc"
    traj = mn.Molecule(mda.Universe(topo, [coords, coords]), name="multi")
    n_frames = traj.universe.trajectory.n_frames
    blend_path = tmp_path / "test.blend"

    session = mn.session.get_session()
    with chdir(tmp_path):
        session.pickle(blend_path)
        session.clear()
        session.load(blend_path)

    restored = session.get(traj.uuid)
    filenames = [Path(f) for f in restored.universe.trajectory.filenames]
    assert len(filenames) == 2
    assert all(f.is_absolute() for f in filenames)
    assert restored.universe.trajectory.n_frames == n_frames


def test_entity_blender_properties(session: MNSession, universe):
    assert len(session.entities) == 0
    t1 = mn.Molecule(universe, name="u1")
    # entity is tracked in the session, keyed by uuid
    assert len(session.entities) == 1
    assert t1.uuid in session.entities
    obj = bpy.data.objects["u1"]
    # entity type is stored on the object properties
    assert obj.mn.entity_type == "molecule"
    assert obj.uuid == t1.uuid
    # visibility is driven through the object properties
    assert obj.mn.visible
    assert obj.visible_get()
    obj.mn.visible = False
    assert not obj.visible_get()
    # removal drops the entity from the session
    session.remove_entity(t1.uuid)
    assert len(session.entities) == 0


def test_remove_tolerates_an_already_deleted_object(fetch):
    """Deleting the object from the outliner leaves the entity behind.

    Removing it must still drop the entity rather than raising
    ``LinkedObjectError`` on the missing object.
    """
    session = mn.session.get_session()
    session.clear()
    mol = fetch("4ozs")
    bpy.data.objects.remove(mol.object, do_unlink=True)

    session.remove(mol.uuid)
    assert mol.uuid not in session.entities


def test_remove_unknown_uuid_raises():
    with pytest.raises(ValueError, match="No entity with UUID"):
        mn.session.get_session().remove("not-a-uuid")


def test_session_saves_when_an_object_was_deleted(fetch, tmp_path):
    """A molecule deleted from the outliner used to abort the session save.

    The .blend still wrote, so the failure was silent - every other entity in
    the file lost its session state on reload.
    """
    session = mn.session.get_session()
    session.clear()
    deleted = fetch("4ozs")
    kept = fetch("1cbs")
    bpy.data.objects.remove(deleted.object, do_unlink=True)

    filepath = tmp_path / "deleted.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(filepath))
    assert session.stashpath(filepath).exists()

    session.clear()
    session.load(filepath)
    assert kept.uuid in session.entities
    assert deleted.uuid not in session.entities


def test_opening_another_file_drops_the_previous_entities(fetch, tmp_path):
    "Entities from a closed file used to linger and raise LinkedObjectError."
    session = mn.session.get_session()
    session.clear()
    fetch("4ozs")
    with_molecule = tmp_path / "with_molecule.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(with_molecule))

    empty = tmp_path / "empty.blend"
    bpy.ops.wm.read_homefile(app_template="")
    bpy.ops.wm.save_as_mainfile(filepath=str(empty))

    bpy.ops.wm.open_mainfile(filepath=str(with_molecule))
    assert session.n_items == 1

    bpy.ops.wm.open_mainfile(filepath=str(empty))
    assert session.n_items == 0


def test_session_load_skips_entity_that_fails_to_restore(tmp_path):
    "One entity failing to restore used to lose every entity in the session"
    import shutil

    session = mn.session.get_session()
    session.clear()
    topo = data_dir / "md_ppr/box.gro"
    coords = tmp_path / "moved.xtc"
    shutil.copy(data_dir / "md_ppr/first_5_frames.xtc", coords)
    good = mn.Molecule.load(topo, data_dir / "md_ppr/first_5_frames.xtc")
    bad = mn.Molecule.load(topo, coords, name="moved")
    blend_path = tmp_path / "test.blend"

    session.pickle(blend_path)
    # release the file handle so the file can be removed on Windows
    bad.universe.trajectory.close()
    coords.unlink()
    session.clear()
    with pytest.warns(UserWarning, match="Could not restore `moved`.*Reload"):
        session.load(blend_path)

    assert session.get(good.uuid) is not None
    assert session.get(bad.uuid) is None


def test_session_loads_legacy_single_pickle(tmp_path):
    "Sessions saved before entities were pickled individually still load"
    import pickle

    session = mn.session.get_session()
    session.clear()
    mol = mn.Molecule.load(data_dir / "1cd3.cif")
    blend_path = tmp_path / "test.blend"
    with open(session.stashpath(blend_path), "wb") as f:
        pickle.dump(session, f)

    session.clear()
    session.load(blend_path)
    assert isinstance(session.get(mol.uuid), mn.Molecule)
