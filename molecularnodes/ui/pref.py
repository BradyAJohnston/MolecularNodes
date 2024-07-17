import bpy

from .. import __package__, template


class MN_OT_Template_Install(bpy.types.Operator):
    bl_idname = "mn.template_install"
    bl_label = "Install Template"
    bl_description = "Install the Molecular Nodes startup template file."

    def execute(self, context):
        template.install()
        self.report({"INFO"}, "Installed Molecular Nodes template.")
        return {"FINISHED"}


class MN_OT_Template_Uninstall(bpy.types.Operator):
    bl_idname = "mn.template_uninstall"
    bl_label = "Uninstall Template"
    bl_description = "Uninstall the Molecular Nodes startup template file."

    @classmethod
    def poll(cls, context):
        return template.is_installed()

    def execute(self, context):
        try:
            template.uninstall()
            self.report({"INFO"}, "Uninstalled Molecular Nodes template.")
        except FileNotFoundError:
            self.report({"WARNING"}, "Template not installed.")

        return {"FINISHED"}


class MolecularNodesPreferences(bpy.types.AddonPreferences):
    bl_idname = __package__

    def draw(self, context):
        layout = self.layout
        layout.label(
            text="Install the Molecular Nodes template file, to start Blender with useful default settings"
        )
        row = layout.row()
        if not template.is_installed():
            text = "Install Template"
        else:
            text = "Reinstall Template"

        row.operator("mn.template_install", text=text)
        row.operator("mn.template_uninstall")


CLASSES = [MN_OT_Template_Install, MN_OT_Template_Uninstall, MolecularNodesPreferences]
