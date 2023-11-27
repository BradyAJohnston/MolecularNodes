import bpy
import pathlib
from .. import pkg
from bpy.types import AddonPreferences

install_instructions = "https://bradyajohnston.github.io/MolecularNodes/installation.html#installing-biotite-mdanalysis"
ADDON_DIR = pathlib.Path(__file__).resolve().parent.parent

bpy.types.Scene.pypi_mirror_provider = bpy.props.StringProperty(
    name = 'pypi_mirror_provider', 
    description = 'PyPI Mirror Provider', 
    options = {'TEXTEDIT_UPDATE','LIBRARY_EDITABLE'}, 
    default = 'Default', 
    subtype = 'NONE', 
    search = pkg.get_pypi_mirror_alias,
    )


def apple_silicon_warning(layout):
    box = layout.box()
    box.alert = True
    box.label(text = "On M1/M2/M3 macOS, manual installation is required.")
    box.operator(
        "wm.url_open", text = "Installation Instructions", icon = 'HELP'
    ).url = install_instructions
    

# Defines the preferences panel for the addon, which shows the buttons for 
# installing and reinstalling the required python packages defined in 'requirements.txt'
class MolecularNodesPreferences(AddonPreferences):
    bl_idname = 'molecularnodes'

    def draw(self, context):
        layout = self.layout
        layout.label(text = "Install the required packages for MolecularNodes.")
        
        col_main = layout.column(heading = '', align = False)
        row_import = col_main.row()
        row_import.prop(bpy.context.scene, 'pypi_mirror_provider',text='Set PyPI Mirror')
        
        pkgs = pkg.get_pkgs()
        for package in pkgs.values():
            row = layout.row()
            col = row.column()
            row = col.row()
            if pkg._is_apple_silicon and package['name'] == "MDAnalysis":
                apple_silicon_warning(row)
                row = col.row()
                row.enabled = False
                row.alert = True
            pkg.button_install_pkg(
                layout = row, 
                name = package.get('name'), 
                version = package.get('version'), 
                desc = package.get('desc')
                )
