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
# This includes all global variables that need to be accessable from all files

# ------------------- EXTERNAL MODULES -------------------
import bpy
import sys, os
from bpy.props import FloatProperty, PointerProperty

# ---------------- GLOBAL ADDON LOGGER -------------------
import logging
MolecularNodesAddonLogger = logging.getLogger('MolecularNodes')

global mn_folder
mn_folder = os.path.dirname(__file__)

# clas for molecular nodes addon
class MolecularNodesAddon:
    name = None
    
    # path to the addon directory
    path = bpy.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    tmp_path = bpy.path.abspath(path + "/tmp/")
    libpath = bpy.path.abspath(path + "/lib/")
    logpath = bpy.path.abspath(path + "/logs/")
    
    # external python dependencies
    external_dependencies = [
        					('biotite', 'biotite', '0.35.0'),
        					('MDAnalysis', 'MDAnalysis', '2.2.0'),
             				]
    external_dependencies_installer = False
    
    # append the add-on's path to Blender's python PATH
    sys.path.insert(0, libpath)
    
    # ADDON CHECKS
    # +++++++++++++++++++++++++++++++++++++++
    # check if the specified module can be found in the "lib" directory
    @classmethod
    def is_installed(cls, module, debug=False):
        import importlib.machinery
        import sys
        if sys.version_info.major >= 3 or (sys.version_info.major == 3 and sys.version_info.minor >= 8):
            from importlib.metadata import version

        # extract info
        module_name, install_name, install_version = module

        # try to find the module in the "lib" directory
        module_spec = (importlib.machinery.PathFinder().find_spec(module_name, [cls.libpath]))
        if module_spec:
            if install_version:

                # check if the installed module version fits
                version_comparison = [ a >= b for a,b in zip(list(map(int, version(install_name).split('.'))), list(map(int, install_version.split('.'))))]
                if all(version_comparison):
                    if debug: MolecularNodesAddonLogger.info(" [#] Found module '%s' v.%s." % (module_name, version(install_name)))
                    return True

                else:
                    if debug: MolecularNodesAddonLogger.info(" [#] Found module '%s' v.%s, but require version %s." % (module_name, version(install_name), install_version))
                    return False
            else:
                if debug: MolecularNodesAddonLogger.info(" [#] Found module '%s' v.%s." % (module_name, version(install_name)))
                return True

        if debug: MolecularNodesAddonLogger.info(" [#] Could not find module '%s'." % module_name)
        return False

    # check if all defined dependencies can be found in the "lib" directory
    @classmethod
    def check_dependecies(cls, debug=False):

        # status
        found_all = True

        # are all modules in the packages list available in the "lib" directory?
        for module in cls.external_dependencies:
            if not cls.is_installed(module, debug):
                found_all = False

        return found_all

    # unload all dependencies
    @classmethod
    def unload_dependecies(cls):

        # if the addon is not in installer mode
        if not MolecularNodesAddon.external_dependecies_installer:

            # are all modules in the packages list available in the "lib" directory?
            for module in cls.external_dependecies:

                # get names
                module_name, install_name, install_version = module

                # unload the module
                del sys.modules[module_name]
                #del module_name