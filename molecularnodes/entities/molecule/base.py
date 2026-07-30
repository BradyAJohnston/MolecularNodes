"""Backwards-compatibility shim for the unified ``Molecule`` entity.

``Molecule`` used to be a distinct, biotite-``AtomArray``-backed entity defined here.
It has been unified with the MDAnalysis-``Universe``-backed entity (formerly
``Trajectory``): a single :class:`~molecularnodes.entities.trajectory.base.Molecule`
now represents both single structures and trajectories, using a ``Universe`` as its
data model. Structure files are parsed with biotite and converted to a ``Universe``
via :func:`~molecularnodes.converters.universe_from_atoms`.

This module re-exports the unified class so the old ``entities.molecule.base.Molecule``
import path keeps working. To parse a structure file into a biotite array directly, use
:func:`molecularnodes.entities.molecule.reader.read_structure`.
"""

from ..trajectory.base import Molecule

__all__ = ["Molecule"]
