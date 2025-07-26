from pathlib import Path
import bpy


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


def viewport_tag_redraw() -> None:
    """Tag all the viewport areas for a redraw"""
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type == "VIEW_3D":
                area.tag_redraw()


def get_viewport_region_from_context(context) -> tuple:
    """Get the 3D viewport region and region data from context"""
    region = context.region
    rv3d = None
    if context.space_data is None:
        return (region, rv3d)
    if not context.space_data.region_quadviews:
        rv3d = context.space_data.region_3d
    else:
        # handle quadview case
        if context.area.type != "VIEW_3D" or context.space_data.type != "VIEW_3D":
            return (region, rv3d)
        i = -1
        for region in context.area.regions:
            if region.type == "WINDOW":
                i += 1
                if context.region == region:
                    break
        else:
            return (region, rv3d)
        rv3d = context.space_data.region_quadviews[i]
    return (region, rv3d)
