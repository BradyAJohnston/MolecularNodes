"""Typed access to the properties Molecular Nodes registers on Blender types.

Properties assigned onto Blender types at register time, such as
``bpy.types.Object.mn = PointerProperty(...)`` in ``ui/addon.py``, are unknown
to the bpy stubs. These functions are the one place that reads them, giving
type checkers the registered type.
"""

from typing import TYPE_CHECKING, Any
import bpy

if TYPE_CHECKING:
    from ..session import MNSession
    from ..ui.props import (
        MolecularNodesObjectProperties,
        MolecularNodesSceneProperties,
        TrajectorySelectionItem,
    )


def object_props(obj: bpy.types.Object) -> "MolecularNodesObjectProperties":
    """The ``Object.mn`` properties of an object."""
    return obj.mn  # ty: ignore[unresolved-attribute]


def object_uuid(obj: bpy.types.Object) -> str:
    """The ``Object.uuid`` linking an object to its entity in the session."""
    return obj.uuid  # ty: ignore[unresolved-attribute]


def object_selections(
    obj: bpy.types.Object,
) -> "bpy.types.bpy_prop_collection_idprop[TrajectorySelectionItem]":
    """The ``Object.mn_trajectory_selections`` of a trajectory object."""
    return obj.mn_trajectory_selections  # ty: ignore[unresolved-attribute]


def object_annotations(
    obj: bpy.types.Object,
) -> "bpy.types.bpy_prop_collection_idprop[Any]":
    """The ``Object.mn_annotations`` of an object.

    The item type is created at runtime from the registered annotation types, so
    its attributes are not known statically.
    """
    return obj.mn_annotations  # ty: ignore[unresolved-attribute]


def scene_props(scene: bpy.types.Scene) -> "MolecularNodesSceneProperties":
    """The ``Scene.mn`` properties of a scene."""
    return scene.mn  # ty: ignore[unresolved-attribute]


def scene_session(scene: bpy.types.Scene) -> "MNSession":
    """The ``Scene.MNSession`` holding the entities of the scene."""
    return scene.MNSession  # ty: ignore[unresolved-attribute]
