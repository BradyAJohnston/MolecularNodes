import bpy
from bpy.types import AddonPreferences

from ..template import template_install, template_uninstall
from .. import __package__


class MN_OT_Template_Install(bpy.types.Operator):
    bl_idname = "mn.template_install"
    bl_label = "Install Template File"
    bl_description = "Install the Molecular Nodes startup template file."

    def execute(self, context):
        template_install()
        self.report({"INFO"}, "Installed Molecular Nodes template.")
        return {"FINISHED"}


class MN_OT_Template_Uninstall(bpy.types.Operator):
    bl_idname = "mn.template_uninstall"
    bl_label = "Uninstall Template"
    bl_description = "Uninstall the Molecular Nodes startup template file."

    def execute(self, context):
        template_uninstall()
        self.report({"INFO"}, "Uninstalled Molecular Nodes template.")
        return {"FINISHED"}


class MolecularNodesPreferences(AddonPreferences):
    bl_idname = __package__

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.operator("mn.template_install")
        row.operator("mn.template_uninstall")


CLASSES = [MN_OT_Template_Install, MN_OT_Template_Uninstall, MolecularNodesPreferences]
