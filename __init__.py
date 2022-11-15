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

packages.verify()

# operator that calls the function to import the structure
class MOL_OT_Import_Protein_RCSB(bpy.types.Operator):
    bl_idname = "mol.import_protein_rcsb"
    bl_label = "import_protein_fetch_pdb"
    bl_description = "Download and open a structure from the Protein Data Bank"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return not False

    def execute(self, context):
        open.import_protein_pdb(
            pdb_code = bpy.context.scene.mol_pdb_code, 
            center_molecule = bpy.context.scene.mol_import_center,
            del_solvent = bpy.context.scene.mol_import_del_solvent
            )
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)


def MOL_PT_panel_ui(layout_function, ): 
    col_main = layout_function.column(heading = '', align = False)
    col_main.alert = False
    col_main.enabled = True
    col_main.active = True
    col_main.use_property_split = False
    col_main.use_property_decorate = False
    col_main.scale_x = 1.0
    col_main.scale_y = 1.0
    col_main.alignment = 'Expand'.upper()
    col_main.label(text = "testing", icon_value = 3)
    row_import = col_main.row()
    row_import.prop(bpy.context.scene, 'mol_pdb_code', text='PDB ID', icon_value=0, emboss=True)
    row_import.operator('mol.import_protein_rcsb', text='Download', icon_value=169, emboss=True, depress=False)
    row_options = col_main.row()
    row_options.prop(bpy.context.scene, 'mol_import_center', text='Centre Structre', icon_value=0, emboss=True)
    row_options.prop(bpy.context.scene, 'mol_import_del_solvent', text='Delete Solvent', icon_value=0, emboss=True)

class MOL_PT_panel(bpy.types.Panel):
    bl_label = 'Molecular Nodes'
    bl_idname = 'MOL_PT_panel'
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = 'scene'
    bl_order = 0
    bl_options = {'HEADER_LAYOUT_EXPAND'}
    bl_ui_units_x=0

    @classmethod
    def poll(cls, context):
        return not (False)

    def draw_header(self, context):
        layout = self.layout

    def draw(self, context):
        layout = self.layout
        layout_function = layout
        MOL_PT_panel_ui(layout_function, )



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

    bpy.utils.register_class(MOL_PT_panel)
    bpy.utils.register_class(MOL_OT_Import_Protein_RCSB)


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

    bpy.utils.unregister_class(MOL_PT_panel)
    bpy.utils.unregister_class(MOL_OT_Import_Protein_RCSB)
