# ##### BEGIN GPL LICENSE BLOCK #####
#
#  Copyright © 2022 Brady Johnston
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ##### END GPL LICENSE BLOCK #####

# MODULE DESCRIPTION:
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# This includes everything that is related to the add-on preferences,
# installation of requirements, etc.

# ------------------ INTERNAL MODULES --------------------
from .globals import *
from . import pkg

# ------------------- EXTERNAL MODULES -------------------
import bpy
import os

# ---------------- GLOBAL ADDON LOGGER -------------------

# ------------- DEFINE ADDON PREFERENCES ----------------
# an operator that installs the python dependencies

class MOL_OT_install_dependencies(bpy.types.Operator):
    bl_idname = "mol.install_dependencies"
    bl_label = "Install Dependencies"
    bl_description = "Install the required python packages to enable import."
    bl_options = {'REGISTER', 'INTERNAL'}
    
    def execute(self, context):
        if not pkg.available():
            import datetime
            
            # generate logfile
            logfile_path = os.path.abspath(MolecularNodesAddon.logpath + 'side-packages-install.log')
            logfile = open(logfile_path, 'a')
            
            logfile.write("-----------------------------------" + '\n')
            logfile.write("Installer Started: " + str(datetime.datetime.now()) + '\n')
            logfile.write("-----------------------------------" + '\n')
            
            pkg.install()
            
            logfile.write("###################################" + '\n')
            logfile.write("Installer finished: " + str(datetime.datetime.now()) + '\n')
            logfile.write("###################################" + '\n')
            
            # close the logfile
            logfile.close()
        
        if pkg.available():
            # bpy.context.preferences.addons['MolecularNodesPref'].preferences.packages_available = True
            self.report(
                {'INFO'}, 
                message='Successfully Installed Required Packages'
                )
        else:
            # bpy.context.preferences.addons['MolecularNodesPref'].preferences.packages_available = False
            self.report(
                {'ERROR'}, 
                message='Failed to install required packages. Please check log file: ' + logfile_path
                )
        
        return {'FINISHED'}


# preferences pane for this Addon in the Blender preferences
class MOL_PT_AddonPreferences(bpy.types.AddonPreferences):
    bl_idname = 'MolecularNodesPref'
    packages_available: bpy.props.BoolProperty(name = 'packages_available', default = False)
    
    def draw(self, context):
        layout = self.layout
        
        if not pkg.available():
            row = layout.row()
            row.alert = True
            row.label(text="Need to install required python packages: 'biotite' and 'MDAnalysis'", icon = 'ERROR')
            row.operator("mol.install_dependencies", icon = "PLUS")
        else:
            row = layout.row()
            row.label(text = 'Biotite and MDAnalysis are installed and available.', icon = "LOCKVIEW_ON")
