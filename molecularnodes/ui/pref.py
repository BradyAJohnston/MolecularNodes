import bpy
from bpy.types import AddonPreferences

from ..utils import template_install

install_instructions = "https://bradyajohnston.github.io/MolecularNodes/installation.html#installing-biotite-mdanalysis"


class MN_OT_Install_Template(bpy.types.Operator):
    bl_idname = "mn.install_template"
    bl_label = "Install Template File"
    bl_description = "Install the Molecular Nodes startup template file."

    def exectute(self, context):
        template_install()


class MolecularNodesPreferences(AddonPreferences):
    bl_idname = __package__

    def draw(self, context):
        layout = self.layout
        layout.label(text="testing")
        layout.operator("mn.install_template", text="Install Template")


CLASSES = [MN_OT_Install_Template, MolecularNodesPreferences]
