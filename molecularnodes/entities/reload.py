"""
Reload (relink) the Python-side session entity for an entity object.

When a ``.blend`` file is opened, the Blender objects for molecular entities are
restored, but the Python objects that back them in the ``MNSession`` (holding the
biotite ``AtomArray``, MDAnalysis ``Universe``, etc.) may be gone. Reloading
reconstructs that Python entity from the source data recorded on the object's
``mn`` properties (file path, PDB code, database) and links it back to the
existing Blender object, without rebuilding its geometry.
"""

import bpy
import MDAnalysis as mda
from ..blender.utils import path_resolve
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
    mn = obj.mn
    # a molecule loaded from topology (+ trajectory) files records those paths;
    # one parsed from a structure file or fetched records filepath/code instead
    if mn.filepath_topology:
        return _reload_trajectory(obj)
    if mn.code:
        file_path = StructureDownloader().download(
            code=mn.code, format="bcif", database=mn.database or "rcsb"
        )
    elif mn.filepath:
        file_path = path_resolve(mn.filepath)
    else:
        raise ValueError("No source file or PDB code recorded for this molecule")
    reader = read_structure(file_path)
    universe = universe_from_atoms(reader.array)
    entity = Molecule(universe, create_object=False)
    entity.object = obj
    return entity


def _reload_trajectory(obj: bpy.types.Object) -> Molecule:
    mn = obj.mn
    path_topo = path_resolve(mn.filepath_topology)
    path_traj = path_resolve(mn.filepath_trajectory)
    if "oxdna" in mn.entity_type:
        universe = mda.Universe(
            path_topo,
            path_traj,
            topology_format=OXDNAParser,
            format=OXDNAReader,
        )
        entity = OXDNA(universe, create_object=False)
    else:
        entity = Molecule.load(path_topo, path_traj, create_object=False)
    entity.object = obj
    entity.set_frame(bpy.context.scene.frame_current)
    return entity


def _reload_density(obj: bpy.types.Object) -> Grids:
    mn = obj.mn
    if not mn.filepath:
        raise ValueError("No source file recorded for this density")
    entity = Grids(file_path=path_resolve(mn.filepath))
    entity.object = obj
    return entity


def _reload_ensemble_star(obj: bpy.types.Object) -> StarFile:
    # StarFile records its source path in the `mn.filepath` property
    return StarFile.from_blender_object(obj)


def _reload_ensemble_cellpack(obj: bpy.types.Object) -> CellPack:
    mn = obj.mn
    if not mn.filepath:
        raise ValueError("No source file recorded for this ensemble")
    entity = CellPack(path_resolve(mn.filepath))
    entity.object = obj
    return entity


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
    """Whether ``obj`` is an entity object that can be reloaded into the session."""
    return obj is not None and obj.mn.entity_type in _RELOADERS


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
    reloader = _RELOADERS.get(obj.mn.entity_type)
    if reloader is None:
        raise ValueError(
            f"Cannot reload an object with entity type '{obj.mn.entity_type}'"
        )
    return reloader(obj)
