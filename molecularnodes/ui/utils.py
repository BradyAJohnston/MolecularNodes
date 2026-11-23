from collections.abc import Callable
from typing import Any
import bpy


def operator_props(
    props: bpy.types.OperatorProperties | None,
) -> bpy.types.OperatorProperties:
    """The properties of an operator button, to fill in.

    ``UILayout.operator()`` is typed as possibly returning None, which it doesn't
    for a registered operator, so wrap its call to set the button's properties.
    """
    assert props is not None
    return props


def menu_draw_func(
    func: Callable[[Any, bpy.types.Context], None],
) -> Callable[..., None]:
    """Mark a function as a draw function for ``Menu.append()`` and ``remove()``.

    Blender calls an appended draw function as ``func(self, context)``, with the
    menu as ``self``, but the bpy stubs type it as ``func(context)``.
    """
    return func


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
        op = operator_props(
            layout.operator("wm.url_open", text="Online Access Docs", icon="URL")
        )
        op.url = "https://docs.blender.org/manual/en/dev/editors/preferences/system.html#bpy-types-preferencessystem-use-online-access"
        layout = layout.column()
        layout.alert = True
        layout.enabled = False

    return layout
