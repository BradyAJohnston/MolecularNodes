import bpy


def check_online_access_for_ui(layout: bpy.types.UILayout) -> bpy.types.UILayout:
    """
    Disable UI without Online Access

    Checks for the online access permissions, and adds a warning and disables following
    UI elements if it fails the check. Returns the UILayout that will have .enabled flag
    set to False, disabling all subsequent uses of the layout.

    Args:
        layout (bpy.types.UILayout): The UILayout element to add the warning and potentially
        disable.

    Returns:
        bpy.types.UILayout: The altered UILayout element, for use in downstream UI
        components.
    """
    if not bpy.app.online_access:
        layout.label(
            text="Online access disabled. Change in Blender's system preferences.",
            icon="ERROR",
        )
        op = layout.operator("wm.url_open", text="Online Access Docs", icon="URL")
        op.url = "https://docs.blender.org/manual/en/dev/editors/preferences/system.html#bpy-types-preferencessystem-use-online-access"
        layout = layout.column()
        layout.alert = True
        layout.enabled = False

    return layout
