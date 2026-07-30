"""Backwards-compatibility shim.

Molecule annotations were unified with the trajectory annotations when both entities
were merged onto a shared MDAnalysis ``Universe`` backend. The annotation classes now
live in :mod:`molecularnodes.entities.trajectory.annotations`; this module re-exports
them so the old ``entities.molecule.annotations`` import path keeps working.
"""

from ..trajectory.annotations import (
    MoleculeAnnotation,
    MoleculeAnnotationManager,
    MoleculeInfo,
)

__all__ = [
    "MoleculeAnnotation",
    "MoleculeAnnotationManager",
    "MoleculeInfo",
]
