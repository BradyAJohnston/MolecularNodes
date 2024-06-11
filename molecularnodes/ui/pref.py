import bpy
import pathlib
from .. import pkg
from ..utils import template_install
from bpy.types import AddonPreferences

install_instructions = "https://bradyajohnston.github.io/MolecularNodes/installation.html#installing-biotite-mdanalysis"
ADDON_DIR = pathlib.Path(__file__).resolve().parent.parent

bpy.types.Scene.pypi_mirror_provider = bpy.props.StringProperty(
    name="pypi_mirror_provider",
    description="PyPI Mirror Provider",
    options={"TEXTEDIT_UPDATE", "LIBRARY_EDITABLE"},
    default="Default",
    subtype="NONE",
    search=pkg.get_pypi_mirror_alias,
)

# Defines the preferences panel for the addon, which shows the buttons for
# installing and reinstalling the required python packages defined in 'requirements.txt'


class MN_OT_Install_Tempalte(bpy.types.Operator):
    bl_idname = "mn.install_template"

    def exectute(self, context):
        template_install()


class MolecularNodesPreferences(AddonPreferences):
    bl_idname = "molecularnodes"

    def draw(self, context):
        layout = self.layout
        layout.operator("mn.install_template", text="Install Template")
        # layout.label(text="Install the required packages for MolecularNodes.")

        # col_main = layout.column(heading="", align=False)
        # row_import = col_main.row()
        # row_import.prop(
        #     bpy.context.scene, "pypi_mirror_provider", text="Set PyPI Mirror"
        # )

        # pkgs = pkg.get_pkgs()
        # for package in pkgs.values():
        #     row = layout.row()
        #     col = row.column()
        #     row = col.row()
        #     pkg.button_install_pkg(
        #         layout=row,
        #         name=package.get("name"),
        #         version=package.get("version"),
        #         desc=package.get("desc"),
        #     )
