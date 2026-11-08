import bpy
import MDAnalysis as mda
import pytest
import molecularnodes as mn
from molecularnodes.blender.registered import object_props, object_uuid
from molecularnodes.session import MNSession
from .constants import data_dir


def test_session_present():
    assert isinstance(mn.session.get_session(), mn.session.MNSession)


def test_persistent_handlers_added():
    load_handlers = [handler.__name__ for handler in bpy.app.handlers.load_post]  # ty: ignore[unresolved-attribute]  # handlers are plain functions
    save_handlers = [handler.__name__ for handler in bpy.app.handlers.save_pre]  # ty: ignore[unresolved-attribute]  # handlers are plain functions
    assert "_save_pre" in save_handlers
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
    assert object_props(obj).filepath != ""
    assert can_reload(obj)

    # simulate a lost session link (e.g. re-opening a .blend without the pickle)
    session.remove_entity(object_uuid(obj))
    assert session.get(object_uuid(obj)) is None

    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.Molecule)
    assert session.get(object_uuid(obj)) is reloaded
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
    assert object_props(obj).filepath_topology != ""
    assert can_reload(obj)

    session.remove_entity(object_uuid(obj))
    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.Molecule)
    assert session.get(object_uuid(obj)) is reloaded


def test_reload_density_from_file():
    from molecularnodes.entities.reload import can_reload, reload_entity

    session = mn.session.get_session()
    density = mn.entities.density.Grids.load(data_dir / "emd_24805.map.gz")
    obj = density.object
    assert object_props(obj).filepath != ""
    assert can_reload(obj)

    session.remove_entity(object_uuid(obj))
    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.entities.density.Density)
    assert session.get(object_uuid(obj)) is reloaded


def test_color_by_styled_density_survives_save_and_reload(fetch, tmp_path):
    "Colours sampled from a styled density's bundled volume after a .blend round trip."
    import numpy as np
    from gridData import Grid
    from molecularnodes.entities.reload import reload_entity

    session = mn.session.get_session()
    session.clear()
    mol = fetch("1cbs")
    positions = mol.named_attribute("position") / mol.world_scale
    origin = np.floor(positions.min(axis=0)) - 10.0
    shape = tuple(np.ceil(np.ptp(positions, axis=0) + 20).astype(int))
    # potential rising along x from -10 to +10 across the molecule
    x = origin[0] + np.arange(shape[0]) - positions[:, 0].mean()
    data = np.broadcast_to((x * 20 / np.ptp(x))[:, None, None], shape).copy()
    Grid(data, origin=origin, delta=1.0).export(str(tmp_path / "ramp.dx"))
    density = mn.entities.density.Grids.load(
        tmp_path / "ramp.dx", style="density_iso_surface"
    )
    mol.add_style("surface", color=density)
    before = mol.named_attribute("Color", evaluate=True)
    # a spread of colours, not a flat result from a missing grid
    assert np.ptp(before[:, 0]) > 0.3

    filepath = tmp_path / "colored.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(filepath))
    session.clear()
    bpy.ops.wm.open_mainfile(filepath=str(filepath))

    restored = session.get(mol.uuid)
    assert isinstance(restored, mn.Molecule)
    assert np.allclose(restored.named_attribute("Color", evaluate=True), before)

    density_obj = bpy.data.objects[density.name]
    session.remove_entity(object_uuid(density_obj))
    reload_entity(density_obj)
    assert np.allclose(restored.named_attribute("Color", evaluate=True), before)


def test_reload_ensemble_star_from_file():
    from molecularnodes.entities.reload import can_reload, reload_entity

    session = mn.session.get_session()
    ensemble = mn.entities.ensemble.StarFile.load(data_dir / "starfile/relion.star")
    obj = ensemble.object
    assert can_reload(obj)

    session.remove_entity(object_uuid(obj))
    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.entities.ensemble.StarFile)
    assert session.get(object_uuid(obj)) is reloaded


def test_reload_ensemble_cellpack_from_file():
    from molecularnodes.entities.reload import can_reload, reload_entity

    session = mn.session.get_session()
    ensemble = mn.entities.ensemble.CellPack.load(
        file_path=data_dir / "cellpack/square1.bcif",
        name="CellPack",
        node_setup=False,
    )
    obj = ensemble.object
    assert object_props(obj).filepath != ""
    assert can_reload(obj)

    session.remove_entity(object_uuid(obj))
    reloaded = reload_entity(obj)
    assert isinstance(reloaded, mn.entities.ensemble.CellPack)
    assert session.get(object_uuid(obj)) is reloaded


def save_and_reopen(blend_path) -> MNSession:
    "Save the file and open it again, rebuilding the session from the objects"
    session = mn.session.get_session()
    bpy.ops.wm.save_as_mainfile(filepath=str(blend_path))
    session.clear()
    bpy.ops.wm.open_mainfile(filepath=str(blend_path))
    return session


def test_session_roundtrip_cellpack(tmp_path):
    ensemble = mn.entities.ensemble.CellPack.load(
        data_dir / "cellpack/square1.bcif", node_setup=False
    )

    session = save_and_reopen(tmp_path / "test.blend")
    restored = session.get(ensemble.uuid)
    assert isinstance(restored, mn.entities.ensemble.CellPack)
    assert restored.name == ensemble.name
    assert restored.instance_collection is not None


@pytest.fixture()
def session():
    return mn.session.get_session()


@pytest.fixture()
def universe():
    topo = data_dir / "md_ppr/box.gro"
    traj = data_dir / "md_ppr/first_5_frames.xtc"
    return mda.Universe(topo, traj)


def test_saving_makes_source_paths_relative(tmp_path):
    "A project folder moved as a whole still restores its entities"
    import shutil

    project = tmp_path / "project"
    (project / "data").mkdir(parents=True)
    for name in ["box.gro", "first_5_frames.xtc"]:
        shutil.copy(data_dir / "md_ppr" / name, project / "data" / name)
    coords = project / "data/first_5_frames.xtc"
    traj = mn.Molecule(
        mda.Universe(project / "data/box.gro", [coords, coords]), name="multi"
    )
    n_frames = traj.universe.trajectory.n_frames
    bpy.ops.wm.save_as_mainfile(filepath=str(project / "test.blend"))

    props = object_props(traj.object)
    assert props.filepath_topology == "//data/box.gro"
    assert props.filepath_trajectory == "//data/first_5_frames.xtc"
    assert props.filepaths_trajectory == ["//data/first_5_frames.xtc"] * 2

    # release the file handles so the folder can be moved on Windows
    traj.universe.trajectory.close()
    moved = tmp_path / "moved"
    project.rename(moved)
    session = mn.session.get_session()
    session.clear()
    bpy.ops.wm.open_mainfile(filepath=str(moved / "test.blend"))

    restored = session.get(traj.uuid)
    assert isinstance(restored, mn.Molecule)
    assert restored.universe.trajectory.n_frames == n_frames


def test_session_roundtrip_multi_file_trajectory(tmp_path):
    from pathlib import Path

    topo = data_dir / "md_ppr/box.gro"
    coords = data_dir / "md_ppr/first_5_frames.xtc"
    traj = mn.Molecule(mda.Universe(topo, [coords, coords]), name="multi")
    n_frames = traj.universe.trajectory.n_frames

    session = save_and_reopen(tmp_path / "test.blend")

    restored = session.get(traj.uuid)
    assert isinstance(restored, mn.Molecule)
    filenames = [Path(f) for f in restored.universe.trajectory.filenames]
    assert len(filenames) == 2
    assert all(f.is_absolute() for f in filenames)
    assert restored.universe.trajectory.n_frames == n_frames


def test_session_roundtrip_keeps_annotations(tmp_path):
    "Annotations used to be lost when an entity was rebuilt from its object"
    mol = mn.Molecule.load(data_dir / "1cd3.cif")
    mol.annotations.add_label_3d(location=(0, 0, 0), text="label")

    session = save_and_reopen(tmp_path / "test.blend")

    restored = session.get(mol.uuid)
    assert restored is not None
    assert len(restored.annotations) == 1


def test_entity_blender_properties(session: MNSession, universe):
    assert len(session.entities) == 0
    t1 = mn.Molecule(universe, name="u1")
    # entity is tracked in the session, keyed by uuid
    assert len(session.entities) == 1
    assert t1.uuid in session.entities
    obj = bpy.data.objects["u1"]
    # entity type is stored on the object properties
    assert object_props(obj).entity_type == "molecule"
    assert object_uuid(obj) == t1.uuid
    # visibility is driven through the object properties
    assert object_props(obj).visible
    assert obj.visible_get()
    object_props(obj).visible = False
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

    save_and_reopen(tmp_path / "deleted.blend")
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

    bpy.ops.wm.save_as_mainfile(filepath=str(blend_path))
    # release the file handle so the file can be removed on Windows
    bad.universe.trajectory.close()
    coords.unlink()
    session.clear()
    with pytest.warns(UserWarning, match="Could not restore `moved`.*Reload"):
        bpy.ops.wm.open_mainfile(filepath=str(blend_path))

    assert session.get(good.uuid) is not None
    assert session.get(bad.uuid) is None


class _RunsOnUnpickle:
    "Records that it was unpickled, as a malicious pickle would run code"

    ran = False

    def __reduce__(self):
        return (_mark_unpickled, ())


def _mark_unpickled():
    _RunsOnUnpickle.ran = True


def test_legacy_session_file_is_not_unpickled(tmp_path):
    "A .MNSession file next to a .blend is never read, but the entities still load"
    import pickle

    mol = mn.Molecule.load(data_dir / "1cd3.cif")
    blend_path = tmp_path / "test.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(blend_path))
    with open(f"{blend_path}.MNSession", "wb") as f:
        pickle.dump(_RunsOnUnpickle(), f)

    session = mn.session.get_session()
    session.clear()
    bpy.ops.wm.open_mainfile(filepath=str(blend_path))

    assert not _RunsOnUnpickle.ran
    assert isinstance(session.get(mol.uuid), mn.Molecule)
