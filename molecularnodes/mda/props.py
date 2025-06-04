"""
Blender properties and associated callbacks
"""

import bpy
from bpy.props import (  # type: ignore
    IntProperty,
    BoolProperty,
    StringProperty,
    PointerProperty,
    CollectionProperty,
)
from . import utils


def _get_universe_visibility(self) -> bool:
    """get callback for universe visibility property"""
    # TODO: Fix objects hidden from the outliner
    return self.get("visible", True)


def _set_universe_visibility(self, visible: bool) -> None:
    """set callback for universe visibility property"""
    self["visible"] = visible
    utils.set_object_visibility(self.object, self.visible)
    # TODO: set visibility of linked components (like density etc later)


def _universe_active_index_callback(self, context: bpy.context) -> None:
    """update callback for universe active_index change"""
    if self.active_index == -1:
        return
    active_object = context.active_object
    universe = context.scene.mn.mda.universes[self.active_index]
    if active_object != universe.object:
        # TODO: Frame this universe into view
        bpy.ops.object.select_all(action="DESELECT")  # deselect all objects
        context.view_layer.objects.active = universe.object  # make active object
        bpy.context.view_layer.update()  # update view layer to reflect changes
        if bpy.context.active_object:
            bpy.context.active_object.select_set(True)  # set as selected object


class mdaUniverseProperties(bpy.types.PropertyGroup):
    """
    Scene level properties for each universe

    Attributes
    ----------
    visible: bool
        bool to determine universe visibility
    object: bpy.types.Object
        pointer to Blender object
    """

    visible: BoolProperty(
        name="visible",
        description="Visibility of the universe",
        default=True,
        get=_get_universe_visibility,
        set=_set_universe_visibility,
    )  # type: ignore
    object: PointerProperty(type=bpy.types.Object)  # type: ignore


class mdaSceneProperties(bpy.types.PropertyGroup):
    """
    Scene level properties for all universes

    Attributes
    ----------
    next_index: int
        unique incrementing index to create key for universes collection
    universes: CollectionProperty
        collection of properties for each universe (scene level)
    active_index: int
        currently selected index in the Universes panel
    """

    next_index: IntProperty(default=0)  # type: ignore
    universes: CollectionProperty(type=mdaUniverseProperties)  # type: ignore
    active_index: IntProperty(
        name="Active universe index",
        default=-1,
        update=_universe_active_index_callback,
    )  # type: ignore


class mdaObjectProperties(bpy.types.PropertyGroup):
    """
    Object level properties for each universe

    Attributes
    ----------
    is_mda_universe: bool
        bool to indicate if object is MDA universe
    universe_key: str
        unique key to lookup universe properties from scene universes collection
    """

    is_mda_universe: bpy.props.BoolProperty(default=False)  # type: ignore
    universe_key: StringProperty(default="")  # type: ignore


CLASSES = [
    mdaUniverseProperties,
    mdaSceneProperties,
    mdaObjectProperties,
]
