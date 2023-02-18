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
    "name"        : "MolecularNodes",
    "author"      : "Brady Johnston", 
    "description" : "Importer and nodes for working with structural biology data in Blender.",
    "blender"     : (3, 4, 0),
    "version"     : (2, 3, 2),
    "location"    : "Scene Properties -> MolecularNodes",
    "warning"     : "",
    "doc_url"     : "https://bradyajohnston.github.io/MolecularNodes/", 
    "tracker_url" : "https://github.com/BradyAJohnston/MolecularNodes/issues", 
    "category"    : "Import"
}

import bpy
from . import pkg
from .pref import *
pkg.verify()
from .load import *
from .ui import *
from .md import *
from .pkg import get_pypi_mirror_alias


def register():
    bpy.types.Scene.pypi_mirror = bpy.props.StringProperty(
        name = 'pypi_mirror', 
        description = 'PyPI Mirror', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = 'Default', 
        subtype = 'NONE', 
        search = get_pypi_mirror_alias,
        )
    bpy.types.Scene.mol_pdb_code = bpy.props.StringProperty(
        name = 'pdb_code', 
        description = 'The 4-character PDB code to download', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = '1bna', 
        subtype = 'NONE', 
        maxlen = 4
        )
    bpy.types.Scene.mol_md_selection = bpy.props.StringProperty(
        name = 'md_selection', 
        description = 'Custom selection string when importing MD simulation. See: "https://docs.mdanalysis.org/stable/documentation_pages/selections.html"', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = 'not (name H* or name OW)', 
        subtype = 'NONE'
        )
    bpy.types.Scene.mol_import_center = bpy.props.BoolProperty(
        name = "mol_import_centre", 
        description = "Move the imported Molecule on the World Origin",
        default = False
        )
    bpy.types.Scene.mol_import_del_solvent = bpy.props.BoolProperty(
        name = "mol_import_del_solvent", 
        description = "Delete the solvent from the structure on import",
        default = True
        )
    bpy.types.Scene.mol_import_include_bonds = bpy.props.BoolProperty(
        name = "mol_import_include_bonds", 
        description = "Include bonds in the imported structure.",
        default = True
        )
    bpy.types.Scene.mol_import_panel_selection = bpy.props.IntProperty(
        name = "mol_import_panel_selection", 
        description = "Import Panel Selection", 
        subtype = 'NONE',
        default = 0
    )
    bpy.types.Scene.mol_import_local_path = bpy.props.StringProperty(
        name = 'path_pdb', 
        description = 'File path of the structure to open', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = '', 
        subtype = 'FILE_PATH', 
        maxlen = 0
        )
    bpy.types.Scene.mol_import_md_topology = bpy.props.StringProperty(
        name = 'path_topology', 
        description = 'File path for the toplogy file for the trajectory', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = '', 
        subtype = 'FILE_PATH', 
        maxlen = 0
        )
    bpy.types.Scene.mol_import_md_trajectory = bpy.props.StringProperty(
        name = 'path_trajectory', 
        description = 'File path for the trajectory file for the trajectory', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = '', 
        subtype = 'FILE_PATH', 
        maxlen = 0
        )
    bpy.types.Scene.mol_import_local_name = bpy.props.StringProperty(
        name = 'mol_name', 
        description = 'Name of the molecule on import', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = '', 
        subtype = 'NONE', 
        maxlen = 0
        )
    bpy.types.Scene.mol_import_md_name = bpy.props.StringProperty(
        name = 'mol_md_name', 
        description = 'Name of the molecule on import', 
        options = {'TEXTEDIT_UPDATE'}, 
        default = '', 
        subtype = 'NONE', 
        maxlen = 0
        )
    bpy.types.Scene.mol_import_md_frame_start = bpy.props.IntProperty(
        name = "mol_import_md_frame_start", 
        description = "Frame start for importing MD trajectory", 
        subtype = 'NONE',
        default = 0
    )
    bpy.types.Scene.mol_import_md_frame_step = bpy.props.IntProperty(
        name = "mol_import_md_frame_step", 
        description = "Frame step for importing MD trajectory", 
        subtype = 'NONE',
        default = 1
    )
    bpy.types.Scene.mol_import_md_frame_end = bpy.props.IntProperty(
        name = "mol_import_md_frame_end", 
        description = "Frame end for importing MD trajectory", 
        subtype = 'NONE',
        default = 50
    )
    bpy.types.Scene.mol_import_default_style = bpy.props.IntProperty(
        name = "mol_import_default_style", 
        description = "Default style for importing molecules.", 
        subtype = 'NONE',
        default = 0
    )

    bpy.utils.register_class(TrajectorySelectionList)
    bpy.utils.register_class(TrajectorySelectionListUI)
    bpy.utils.register_class(TrajectorySelection_OT_NewItem)
    bpy.utils.register_class(TrajectorySelection_OT_DeleteIem)
    
    bpy.types.Scene.trajectory_selection_list = bpy.props.CollectionProperty(
        type = TrajectorySelectionList
    )
    bpy.types.Scene.list_index = bpy.props.IntProperty(
        name = "Index for trajectory selection list.", 
        default = 0
    )
    
    bpy.types.NODE_MT_add.append(mol_add_node_menu)

    bpy.utils.register_class(MOL_PT_panel)
    bpy.utils.register_class(MOL_PT_AddonPreferences)
    bpy.utils.register_class(MOL_MT_Add_Node_Menu)
    bpy.utils.register_class(MOL_MT_Add_Node_Menu_Properties)
    bpy.utils.register_class(MOL_MT_Add_Node_Menu_Styling)
    bpy.utils.register_class(MOL_MT_Add_Node_Menu_Color)
    bpy.utils.register_class(MOL_MT_Add_Node_Menu_Bonds)
    bpy.utils.register_class(MOL_MT_Add_Node_Menu_Selections)
    bpy.utils.register_class(MOL_MT_Add_Node_Menu_Membranes)
    bpy.utils.register_class(MOL_MT_Add_Node_Menu_DNA)
    bpy.utils.register_class(MOL_MT_Add_Node_Menu_Animation)
    bpy.utils.register_class(MOL_MT_Add_Node_Menu_Utilities)
    bpy.utils.register_class(MOL_MT_Add_Node_Menu_Assembly)
    bpy.utils.register_class(MOL_MT_Default_Style)

    bpy.utils.register_class(MOL_OT_Style_Surface_Custom)

    bpy.utils.register_class(MOL_OT_Import_Protein_RCSB)

    bpy.utils.register_class(MOL_OT_Import_Method_Selection)
    bpy.utils.register_class(MOL_OT_Import_Protein_Local)
    bpy.utils.register_class(MOL_OT_Import_Protein_MD)
    bpy.utils.register_class(MOL_OT_Assembly_Bio)
    bpy.utils.register_class(MOL_OT_Default_Style)
    bpy.utils.register_class(MOL_OT_Color_Chain)
    bpy.utils.register_class(MOL_OT_Chain_Selection_Custom)
    bpy.utils.register_class(MOL_OT_Ligand_Selection_Custom)
    
    bpy.utils.register_class(MOL_OT_install_dependencies)
    bpy.utils.register_class(MOL_OT_Add_Custom_Node_Group)

    bpy.utils.register_class(MOL_OT_Residues_Selection_Custom)


def unregister():
    del bpy.types.Scene.mol_pdb_code
    del bpy.types.Scene.mol_md_selection
    del bpy.types.Scene.mol_import_center
    del bpy.types.Scene.mol_import_del_solvent
    del bpy.types.Scene.mol_import_include_bonds
    del bpy.types.Scene.mol_import_panel_selection
    del bpy.types.Scene.mol_import_local_path
    del bpy.types.Scene.mol_import_md_topology
    del bpy.types.Scene.mol_import_md_trajectory
    del bpy.types.Scene.mol_import_local_name
    del bpy.types.Scene.mol_import_md_name
    del bpy.types.Scene.mol_import_md_frame_start
    del bpy.types.Scene.mol_import_md_frame_step
    del bpy.types.Scene.mol_import_md_frame_end
    del bpy.types.Scene.mol_import_default_style
    
    del bpy.types.Scene.trajectory_selection_list
    del bpy.types.Scene.list_index
    
    bpy.types.NODE_MT_add.remove(mol_add_node_menu)
    
    bpy.utils.unregister_class(TrajectorySelectionList)
    bpy.utils.unregister_class(TrajectorySelectionListUI)
    bpy.utils.unregister_class(TrajectorySelection_OT_NewItem)
    bpy.utils.unregister_class(TrajectorySelection_OT_DeleteIem)

    bpy.utils.unregister_class(MOL_PT_panel)
    bpy.utils.unregister_class(MOL_PT_AddonPreferences)
    bpy.utils.unregister_class(MOL_MT_Add_Node_Menu)
    bpy.utils.unregister_class(MOL_MT_Add_Node_Menu_Properties)
    bpy.utils.unregister_class(MOL_MT_Add_Node_Menu_Styling)
    bpy.utils.unregister_class(MOL_MT_Add_Node_Menu_Color)
    bpy.utils.unregister_class(MOL_MT_Add_Node_Menu_Bonds)
    bpy.utils.unregister_class(MOL_MT_Add_Node_Menu_Selections)
    bpy.utils.unregister_class(MOL_MT_Add_Node_Menu_Membranes)
    bpy.utils.unregister_class(MOL_MT_Add_Node_Menu_DNA)
    bpy.utils.unregister_class(MOL_MT_Add_Node_Menu_Animation)
    bpy.utils.unregister_class(MOL_MT_Add_Node_Menu_Utilities)
    bpy.utils.unregister_class(MOL_MT_Add_Node_Menu_Assembly)
    bpy.utils.unregister_class(MOL_MT_Default_Style)
    bpy.utils.unregister_class(MOL_OT_Style_Surface_Custom)
    bpy.utils.unregister_class(MOL_OT_Ligand_Selection_Custom)
    
    bpy.utils.unregister_class(MOL_OT_Import_Protein_RCSB)
    bpy.utils.unregister_class(MOL_OT_Import_Method_Selection)
    bpy.utils.unregister_class(MOL_OT_Import_Protein_Local)
    bpy.utils.unregister_class(MOL_OT_Import_Protein_MD)
    bpy.utils.unregister_class(MOL_OT_Assembly_Bio)
    bpy.utils.unregister_class(MOL_OT_Default_Style)
    bpy.utils.unregister_class(MOL_OT_Color_Chain)
    bpy.utils.unregister_class(MOL_OT_Chain_Selection_Custom)
    
    bpy.utils.unregister_class(MOL_OT_install_dependencies)
    bpy.utils.unregister_class(MOL_OT_Add_Custom_Node_Group)

    bpy.utils.unregister_class(MOL_OT_Residues_Selection_Custom)

if __name__=="__main__":
    register()
