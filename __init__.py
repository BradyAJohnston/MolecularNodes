# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


bl_info = {
    "name"        : "Molecular Nodes2 boogaloo",
    "author"      : "Brady Johnston", 
    "description" : "Some more nodes",
    "blender"     : (3, 3, 0),
    "version"     : (0, 9, 1),
    "location"    : "Perth, Australia",
    "warning"     : "",
    "doc_url"     : "https://bradyajohnston.github.io/MolecularNodes/", 
    "tracker_url" : "https://github.com/BradyAJohnston/MolecularNodes/issues", 
    "category"    : "Molecular"
}

addon_keymaps = {}
_icons = None

import bpy
from .src import open
from .src import packages
from .src.panel import *

packages.verify()

def register():
    global _icons
    _icons = bpy.utils.previews.new()
    bpy.types.Scene.sna_atomium_available = bpy.props.BoolProperty(name='atomium_available', description='', default=False)
    bpy.types.Scene.mol_pdb_code = bpy.props.StringProperty(
        name = 'pdb_code', 
        description = 'The 4-character PDB code to download.', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = '1bna', 
        subtype = 'NONE', 
        maxlen = 4
        )
    bpy.types.Scene.mol_import_center = bpy.props.BoolProperty(
        name = "mol_import_centre", 
        description = "Move the imported Molecule on the World Origin",
        default = True
        )
    bpy.types.Scene.mol_import_del_solvent = bpy.props.BoolProperty(
        name = "mol_import_del_solvent", 
        description = "Delete the solvent from the structure on import",
        default = True
        )
    bpy.types.Scene.mol_import_panel_selection = bpy.props.IntProperty(
        name = "mol_import_panel_selection", 
        description = "Import Panel Selection", 
        subtype = 'NONE',
        default = 0
    )
    bpy.types.Scene.mol_import_local_path = bpy.props.StringProperty(
        name = 'pdb_path', 
        description = 'File path of the structure to open', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = '', 
        subtype = 'FILE_PATH', 
        maxlen = 0
        )
    
    bpy.types.Scene.mol_import_md_topology = bpy.props.StringProperty(
        name = 'pdb_path', 
        description = 'File path for the toplogy file for the trajectory', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = '', 
        subtype = 'FILE_PATH', 
        maxlen = 0
        )
    
    bpy.types.Scene.mol_import_md_trajectory = bpy.props.StringProperty(
        name = 'pdb_path', 
        description = 'File path for the trajectory file for the trajectory', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = '', 
        subtype = 'FILE_PATH', 
        maxlen = 0
        )


    bpy.utils.register_class(MOL_PT_panel)
    bpy.utils.register_class(MOL_OT_Import_Protein_RCSB)
    bpy.utils.register_class(MOL_OT_Import_Method_Selection)


def unregister():
    global _icons
    bpy.utils.previews.remove(_icons)
    wm = bpy.context.window_manager
    kc = wm.keyconfigs.addon
    for km, kmi in addon_keymaps.values():
        km.keymap_items.remove(kmi)
        addon_keymaps.clear()
    
    
    del bpy.types.Scene.mol_pdb_code
    del bpy.types.Scene.mol_import_center
    del bpy.types.Scene.mol_import_del_solvent
    del bpy.types.Scene.mol_import_panel_selection
    del bpy.types.Scene.mol_import_md_topology
    del bpy.types.Scene.mol_import_md_trajectory

    bpy.utils.unregister_class(MOL_PT_panel)
    bpy.utils.unregister_class(MOL_OT_Import_Protein_RCSB)
    bpy.utils.unregister_class(MOL_OT_Import_Method_Selection)
