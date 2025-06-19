from pathlib import Path
import bpy
from bpy.app.handlers import persistent  # type: ignore


def path_resolve(path: str | Path) -> Path:
    if isinstance(path, str):
        return Path(bpy.path.abspath(path))
    elif isinstance(path, Path):
        return Path(bpy.path.abspath(str(path)))
    else:
        raise ValueError(f"Unable to resolve path: {path}")


def set_object_visibility(object: bpy.types.Object, visible: bool) -> None:
    """Set visibility of Blender object"""
    if object.name not in bpy.context.view_layer.objects:
        return
    try:
        object.hide_set(not visible)
        # obj.hide_viewport = hide
        # icon incorrect, but causes blender bug for volumes with above
        object.hide_render = not visible
    except RuntimeError:
        # Keyframing object visibility can at times lead to:
        # RuntimeError: Object can't be hidden because it is not in View Layer 'ViewLayer'!
        pass


def _msgbus_active_object_callback(*args):
    """Callback when active object changes"""
    scene = bpy.context.scene
    object = bpy.context.active_object
    # object can be None if it active but hidden
    if object and object.uuid in scene.MNSession.entities:
        scene.mn.entities_active_index = bpy.data.objects.find(object.name)


def _subscribe_to_active_object_changes():
    """Subscribe to active object changes"""
    subscribe_to = (bpy.types.LayerObjects, "active")
    # all our msgbus subscriptions will use MNSession as owner
    # so that they can easily be unsubscribed during unregister()
    bpy.msgbus.subscribe_rna(
        key=subscribe_to,
        owner=bpy.context.scene.MNSession,
        args=(),
        notify=_msgbus_active_object_callback,
    )


@persistent
def _subscribe_to_msgbus(dummy):
    """Load handler to re-subscribe to msgbus"""
    # Add all required msgbus subscriptions go here
    _subscribe_to_active_object_changes()  # subscribe to active object changes


def register_msgbus_subscriptions():
    """Subscribe to all required msgbus updates"""
    # This method is called during register()
    _subscribe_to_msgbus(None)
    # create a load handler to re-subscribe on file loads
    bpy.app.handlers.load_post.append(_subscribe_to_msgbus)


def unregister_msgbus_subscriptions():
    """Clear all msgbus subscriptions"""
    # This method is called during unregister()
    bpy.msgbus.clear_by_owner(bpy.context.scene.MNSession)
    bpy.app.handlers.load_post.remove(_subscribe_to_msgbus)
