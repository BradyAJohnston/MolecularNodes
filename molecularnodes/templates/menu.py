from pathlib import Path
import bpy

_folder_map = {
    "API": "api",
    "Annotations": "annotations",
}


class MN_MT_templates(bpy.types.Menu):
    bl_idname = "MN_MT_templates"
    bl_label = "Molecular Nodes"

    def draw(self, context):
        layout: bpy.types.UILayout = self.layout
        assert layout is not None

        for i, (label, folder) in enumerate(_folder_map.items()):
            path = Path(__file__).parent / folder
            if not path.is_dir():
                continue
            if i > 0:
                layout.separator()
            layout.label(text=label)
            self.path_menu(
                searchpaths=[path],
                operator="text.open",
                props_default={"internal": True},
            )


def _draw_templates_menu(self, context):
    """Draw method for the main menu item"""
    self.layout.menu("MN_MT_templates")


def register_templates_menu():
    bpy.utils.register_class(MN_MT_templates)
    bpy.types.TEXT_MT_templates.append(_draw_templates_menu)


def unregister_templates_menu():
    bpy.types.TEXT_MT_templates.remove(_draw_templates_menu)
    bpy.utils.unregister_class(MN_MT_templates)
