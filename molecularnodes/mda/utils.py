"""
Utils
"""

import bpy
from bpy.app.handlers import persistent  # type: ignore


def is_mda_universe(object: bpy.types.Object) -> bool:
    """Check if blender object is MDA universe"""
    return object and object.mn.mda and object.mn.mda.is_mda_universe


def _msgbus_active_object_callback(*args):
    """Callback when active object changes"""
    sprop = bpy.context.scene.mn.mda
    active_object = bpy.context.active_object
    index = -1
    if is_mda_universe(active_object):
        # find index of universe in scene universes collection
        index = sprop.universes.find(active_object.mn.mda.universe_key)
    sprop.active_index = index  # updates selection in Universes panel


def _subscribe_to_active_object():
    """Actual subscribe to active object changes"""
    subscribe_to = (bpy.types.LayerObjects, "active")
    bpy.msgbus.subscribe_rna(
        key=subscribe_to,
        owner=bpy,
        args=(),
        notify=_msgbus_active_object_callback,
    )


@persistent
def active_object_load_handler(dummy):
    """Load handler to re-subscribe to active object changes"""
    _subscribe_to_active_object()


def subscribe_to_active_object():
    """Subscribe to active object changes"""
    _subscribe_to_active_object()  # subscribe
    # create a load handler to re-subscribe on file loads
    bpy.app.handlers.load_post.append(active_object_load_handler)


def set_object_visibility(object: bpy.types.Object, visible: bool) -> None:
    """Set visibility of Blender object"""
    if object.name not in bpy.context.view_layer.objects:
        return
    try:
        object.hide_set(not visible)
        # obj.hide_viewport = hide
        # icon incorrect, but causes blender bug for volumes with above
        object.hide_render = not visible
    except:
        # Keyframing object visibility can at times lead to:
        # RuntimeError: Object can't be hidden because it is not in View Layer 'ViewLayer'!
        pass
