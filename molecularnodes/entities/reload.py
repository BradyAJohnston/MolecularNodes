"""
Reload (relink) the Python-side session entity for an entity object.

When a ``.blend`` file is opened, the Blender objects for molecular entities are
restored, but the Python objects that back them in the ``MNSession`` (holding the
biotite ``AtomArray``, MDAnalysis ``Universe``, etc.) may be gone. Reloading
reconstructs that Python entity from the source data recorded on the object's
``mn`` properties (file path, PDB code, database) and links it back to the
existing Blender object, without rebuilding its geometry.
"""

from io import BytesIO
from pathlib import Path
import bpy
import MDAnalysis as mda
from ..blender.registered import object_props, object_uuid, scene_session
from ..blender.utils import path_resolve, resolve_file_path
from ..converters import universe_from_atoms
from ..download import StructureDownloader
from .base import EntityType, MolecularEntity
from .density.grids import Grids
from .ensemble.cellpack import CellPack
from .ensemble.star import StarFile
from .molecule import OXDNA, Molecule
from .molecule.oxdna import OXDNAParser, OXDNAReader
from .molecule.reader import read_structure


def _reload_molecule(obj: bpy.types.Object) -> Molecule:
    mn = object_props(obj)
    # a molecule loaded from topology (+ trajectory) files records those paths;
    # one parsed from a structure file or fetched records filepath/code instead
    if mn.filepath_topology:
        return _reload_trajectory(obj)
    if mn.code:
        file_path = StructureDownloader().download(
            code=mn.code, format="bcif", database=mn.database or "rcsb"
        )
        # binary CIF is never downloaded as text
        assert isinstance(file_path, (Path, BytesIO))
    elif mn.filepath:
        file_path = path_resolve(mn.filepath)
    else:
        raise ValueError("No source file or PDB code recorded for this molecule")
    reader = read_structure(file_path, data_block=mn.data_block or None)
    universe = universe_from_atoms(reader.array, frame_data=reader.model_data())
    entity = Molecule(universe, create_object=False)
    entity.object = obj
    return entity


def _reload_trajectory(obj: bpy.types.Object) -> Molecule:
    mn = object_props(obj)
    path_topo = resolve_file_path(mn.filepath_topology, "Topology file")
    if mn.entity_type == EntityType.MD_OXDNA:
        # an oxDNA topology carries no positions, so the trajectory is required
        path_traj = resolve_file_path(mn.filepath_trajectory, "Trajectory file")
        universe = mda.Universe(
            path_topo,
            path_traj,
            topology_format=OXDNAParser,
            format=OXDNAReader,
        )
        entity = OXDNA(universe, create_object=False)
    elif mn.filepath_trajectory:
        # a multi-file trajectory records every coordinate file, not just the first
        paths_traj = [
            resolve_file_path(path, "Trajectory file")
            for path in mn.filepaths_trajectory or [mn.filepath_trajectory]
        ]
        entity = Molecule.load(
            path_topo,
            paths_traj if len(paths_traj) > 1 else paths_traj[0],
            create_object=False,
        )
    else:
        # a topology-only universe (e.g. a `.gro` with its own coordinates)
        entity = Molecule(mda.Universe(path_topo), create_object=False)
    entity.object = obj
    scene = bpy.context.scene
    assert scene is not None
    entity.set_frame(scene.frame_current)
    return entity


def _reload_density(obj: bpy.types.Object) -> Grids:
    mn = object_props(obj)
    if not mn.filepath:
        raise ValueError("No source file recorded for this density")
    entity = Grids(file_path=path_resolve(mn.filepath))
    entity.object = obj
    return entity


def _reload_ensemble_star(obj: bpy.types.Object) -> StarFile:
    # StarFile records its source path in the `mn.filepath` property
    mn = object_props(obj)
    if not mn.filepath:
        raise ValueError("No source file recorded for this ensemble")
    entity = StarFile(path_resolve(mn.filepath))
    entity.object = obj
    return entity


def _reload_ensemble_cellpack(obj: bpy.types.Object) -> CellPack:
    mn = object_props(obj)
    if not mn.filepath:
        raise ValueError("No source file recorded for this ensemble")
    entity = CellPack(path_resolve(mn.filepath))
    entity.object = obj
    collection = _instance_collection(obj)
    if collection is not None:
        entity.instance_collection = collection
    return entity


def _instance_collection(obj: bpy.types.Object) -> bpy.types.Collection | None:
    """The collection of instances an ensemble object's node tree instances.

    The ensemble only remembers the collection by name in Python, but its tree
    references the collection directly, which survives the collection being renamed.
    """
    for mod in obj.modifiers:
        if not isinstance(mod, bpy.types.NodesModifier) or mod.node_group is None:
            continue
        for node in mod.node_group.nodes:
            for socket in node.inputs:
                if (
                    isinstance(socket, bpy.types.NodeSocketCollection)
                    and socket.default_value is not None
                ):
                    return socket.default_value
    return None


# streaming trajectories are backed by a live IMD connection with no source file,
# so they are intentionally absent and cannot be reloaded
_RELOADERS = {
    EntityType.MOLECULE.value: _reload_molecule,
    EntityType.MD_OXDNA.value: _reload_trajectory,
    EntityType.DENSITY.value: _reload_density,
    EntityType.ENSEMBLE_STAR.value: _reload_ensemble_star,
    EntityType.ENSEMBLE_CELLPACK.value: _reload_ensemble_cellpack,
}


def can_reload(obj: bpy.types.Object) -> bool:
    """Whether ``obj`` is an entity object that can be reloaded into the session.

    The object has to record a source to reload from: the objects a CellPack
    instances share its entity type, but have no source of their own.
    """
    if obj is None:
        return False
    mn = object_props(obj)
    return mn.entity_type in _RELOADERS and bool(
        mn.filepath or mn.code or mn.filepath_topology
    )


def reload_entity(obj: bpy.types.Object) -> MolecularEntity:
    """
    Reconstruct the session entity for an entity object and link it to ``obj``.

    Parameters
    ----------
    obj : bpy.types.Object
        The entity object to reload, using the source recorded on ``obj.mn``.

    Returns
    -------
    MolecularEntity
        The reloaded entity, now registered in the session and linked to ``obj``.

    Raises
    ------
    ValueError
        If the object's entity type cannot be reloaded, or the source data
        needed to reload it is missing.
    """
    entity_type = object_props(obj).entity_type
    reloader = _RELOADERS.get(entity_type)
    if reloader is None:
        raise ValueError(f"Cannot reload an object with entity type '{entity_type}'")
    entity = reloader(obj)
    # the annotation manager is created before the entity is linked to its
    # object, so it couldn't read the annotations stored on the object then
    annotations = getattr(entity, "annotations", None)
    if annotations is not None and len(annotations) == 0:
        annotations._restore_annotation_instances_from_props()
    return entity


def restore_entity(obj: bpy.types.Object) -> MolecularEntity:
    """
    Rebuild the session entity of ``obj`` after a file was opened, keeping its uuid.

    Unlike [](`~reload_entity`), which gives the object a new uuid so that an entity
    still held from before can't take the object back, the entity keeps the uuid
    the object was saved with. Nothing from before the file was opened can still
    hold it, and anything that looks the entity up by that uuid - a script, or the
    session - keeps finding it.

    Parameters
    ----------
    obj : bpy.types.Object
        The entity object to restore, using the source recorded on ``obj.mn``.

    Returns
    -------
    MolecularEntity
        The restored entity, registered in the session under the object's uuid.
    """
    uuid = object_uuid(obj)
    entity = reload_entity(obj)
    if uuid and entity.uuid != uuid:
        scene = bpy.context.scene
        assert scene is not None
        session = scene_session(scene)
        session.entities.pop(entity.uuid, None)
        entity._uuid = uuid
        entity.object = obj
        # every reloader returns a Molecule, Ensemble or Density
        session.register_entity(entity)  # ty: ignore[invalid-argument-type]
    return entity
