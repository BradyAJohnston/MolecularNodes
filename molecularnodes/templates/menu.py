from pathlib import Path
import bpy

# List any new subdirs here
_subdirs = ["api", "annotations"]


def _get_submenu_idname(name):
    """Return unique bl_idname for a submenu item"""
    return f"MN_MT_templates_{name}"


def _create_dynamic_menu_class(subdir):
    """Create a dynamic menu class"""
    path = Path(__file__).parent / subdir

    def draw(self, context):
        self.path_menu(
            searchpaths=[path], operator="text.open", props_default={"internal": True}
        )

    attributes = dict(bl_idname=_get_submenu_idname(subdir), bl_label=subdir, draw=draw)
    return type(subdir, (bpy.types.Menu,), attributes)


class MN_MT_templates(bpy.types.Menu):
    bl_idname = "MN_MT_templates"
    bl_label = "Molecular Nodes"

    def draw(self, context):
        for subdir in _subdirs:
            self.layout.menu(_get_submenu_idname(subdir), text=subdir.title())


def _draw_templates_menu(self, context):
    """Draw method for the main menu item"""
    self.layout.menu("MN_MT_templates")


def register_templates_menu():
    # register submenu classes
    for cls in _submenu_classes:
        bpy.utils.register_class(cls)
    # register main menu class
    bpy.utils.register_class(MN_MT_templates)
    # add to Text Editor > Templates menu
    bpy.types.TEXT_MT_templates.append(_draw_templates_menu)


def unregister_templates_menu():
    # remove from Text Editor > Templates menu
    bpy.types.TEXT_MT_templates.remove(_draw_templates_menu)
    # unregister main menu class
    bpy.utils.unregister_class(MN_MT_templates)
    # unregister submenu classes
    for cls in reversed(_submenu_classes):
        bpy.utils.unregister_class(cls)


_submenu_classes = []
for subdir in _subdirs:
    _submenu_classes.append(_create_dynamic_menu_class(subdir))
