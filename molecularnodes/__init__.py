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

from .io.md import MN_OT_Import_Protein_MD
from .io.star import MN_OT_Import_Star_File
from .io.wwpdb import MN_OT_Import_wwPDB
from .io.local import MN_OT_Import_Protein_Local
from .io.dna import MN_OT_Import_OxDNA_Trajectory
from .io.density import MN_OT_Import_Map
from .io.cellpack import MN_OT_Import_Cell_Pack
import bpy

from .utils import template_install
from .props import MolecularNodesObjectProperties
from .ui.node_menu import MN_add_node_menu
from .ui.panel import MN_PT_panel
from .io.parse.mda import _rejuvenate_universe, _sync_universe
from .io.parse.star import _rehydrate_ensembles
# from .io import ops_io
# from .ui import ops_ui

from .ui.ops import (
    MN_OT_Add_Custom_Node_Group,
    MN_OT_Color_Custom,
    MN_OT_selection_custom,
    MN_OT_Residues_Selection_Custom,
    MN_OT_Change_Style
)

ops_ui = [
    MN_OT_Add_Custom_Node_Group,
    MN_OT_Color_Custom,
    MN_OT_selection_custom,
    MN_OT_Residues_Selection_Custom,
    MN_OT_Change_Style
]


ops_io = [
    MN_OT_Import_Cell_Pack,
    MN_OT_Import_Map,
    MN_OT_Import_OxDNA_Trajectory,
    MN_OT_Import_Protein_Local,
    MN_OT_Import_Protein_MD,
    MN_OT_Import_Star_File,
    MN_OT_Import_wwPDB
]

all_classes = ops_ui + ops_io

universe_funcs = [_sync_universe, _rejuvenate_universe]


def register():
    bpy.utils.register_class(MN_PT_panel)
    bpy.types.NODE_MT_add.append(MN_add_node_menu)

    # register all of the import operators
    for op in all_classes:
        bpy.utils.register_class(op)
    bpy.utils.register_class(MolecularNodesObjectProperties)
    bpy.types.Object.mn = bpy.props.PointerProperty(
        type=MolecularNodesObjectProperties
    )
    for func in universe_funcs:
        try:
            bpy.app.handlers.load_post.append(func)
        except ValueError as e:
            print(f"Filaed to append {func}, error: {e}.")
    template_install()


def unregister():
    bpy.utils.unregister_class(MN_PT_panel)
    bpy.types.NODE_MT_add.remove(MN_add_node_menu)

    # unregister all of the import operator classes
    for op in all_classes:
        bpy.utils.unregister_class(op)

    try:
        bpy.utils.unregister_class(MolecularNodesObjectProperties)
        del bpy.types.Object.mn
        for func in universe_funcs:
            try:
                bpy.app.handlers.load_post.remove(func)
            except ValueError as e:
                print(f"Failed to remove {func}, error: {e}.")
    except RuntimeError:
        pass


# can't register the add-on when these are uncommnted, but they do fix the issue
# of having to call register() when running a script
# unregister()
# register()

# # register won't be called when MN is run as a module
bpy.app.handlers.load_post.append(_rejuvenate_universe)
bpy.app.handlers.load_post.append(_rehydrate_ensembles)
bpy.app.handlers.save_post.append(_sync_universe)
