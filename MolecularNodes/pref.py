import bpy
from . import pkg
from bpy.types import AddonPreferences
from . import ui

# Defines the preferences panel for the addon, which shows the buttons for 
# installing and reinstalling the required python packages defined in 'requirements.txt'
class MolecularNodesPreferences(AddonPreferences):
    bl_idname = 'MolecularNodes'

    def draw(self, context):
        layout = self.layout
        layout.label(text = "Install the required packages for MolecularNodes.")
        pkgs = pkg.get_pkgs()
        for package in pkgs.values():
            row = layout.row()
            ui.button_install_pkg(layout = row, package = package.get('name'))
            if pkg.is_apple_silicon() and package.get('name') == "MDAnalysis":
                row.enabled = False
                if not pkg.is_current('MDAnalysis'):
                    row.enabled = False
                    box = layout.box()
                    box.alert = True
                    box.label(text = "On M1/M2 macOS machines, extra install steps are required.")
                    box.operator(
                        "wm.url_open", text = "Installation Instructions", icon = 'HELP'
                    ).url = "https://bradyajohnston.github.io/MolecularNodes/installation.html#installing-biotite-mdanalysis"