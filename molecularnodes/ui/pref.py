from typing import Annotated, cast
import bpy
from bpy.props import BoolProperty, StringProperty
from .. import __package__ as _addon_package
from .. import assets, download
from ..blender.utils import annotated_props

assert _addon_package is not None
_ADDON_KEY: str = _addon_package


def addon_preferences(
    context: bpy.types.Context | None = None,
) -> "MolecularNodesPreferences | None":
    if context is None:
        context = bpy.context
    preferences = context.preferences
    if preferences is None:
        return None
    try:
        addon = preferences.addons[_ADDON_KEY]
    except KeyError:
        key = "bl_ext.vscode_development.molecularnodes"
        if key not in preferences.addons:
            return None
        addon = preferences.addons[key]
    return cast("MolecularNodesPreferences | None", addon.preferences)


class MN_OT_Template_Install(bpy.types.Operator):
    bl_idname = "mn.template_install"
    bl_label = "Install Template"
    bl_description = "Install the Molecular Nodes startup template file."

    def execute(self, context):
        assets.template.install()
        self.report({"INFO"}, "Installed Molecular Nodes template.")
        return {"FINISHED"}


class MN_OT_Template_Uninstall(bpy.types.Operator):
    bl_idname = "mn.template_uninstall"
    bl_label = "Uninstall Template"
    bl_description = "Uninstall the Molecular Nodes startup template file."

    @classmethod
    def poll(cls, context):
        return assets.template.is_installed()

    def execute(self, context):
        try:
            assets.template.uninstall()
            self.report({"INFO"}, "Uninstalled Molecular Nodes template.")
        except FileNotFoundError:
            self.report({"WARNING"}, "Template not installed.")

        return {"FINISHED"}


@annotated_props
class MolecularNodesPreferences(bpy.types.AddonPreferences):
    bl_idname = _ADDON_KEY

    cache_dir: Annotated[
        str,
        StringProperty(
            name="Cache Directory",
            description="Where to store the structures downloaded from the Protein Data Bank",
            default=str(download.CACHE_DIR),
            subtype="DIR_PATH",
        ),
    ]

    cache_download: Annotated[
        bool,
        BoolProperty(
            name="Cache Downloads",
            description="Cache downloaded files from the Protein Data Bank",
            default=True,
        ),
    ]

    def draw(self, context):
        layout = self.layout
        assert layout is not None
        layout.label(
            text="Where and if to store downloaded files for faster subsequent loading:"
        )
        row = layout.row()
        row.prop(self, "cache_download", text="")
        col = row.column()
        col.prop(self, "cache_dir")
        col.enabled = self.cache_download
        layout.label(
            text="Install the Molecular Nodes template file, to start Blender with useful default settings"
        )
        row = layout.row()
        if not assets.template.is_installed():
            text = "Install Template"
        else:
            text = "Reinstall Template"

        row.operator("mn.template_install", text=text)
        row.operator("mn.template_uninstall")


CLASSES = [MN_OT_Template_Install, MN_OT_Template_Uninstall, MolecularNodesPreferences]
