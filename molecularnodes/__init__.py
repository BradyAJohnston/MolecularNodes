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

import bpy

from .io import ops_io
from .io.md import TrajectorySelectionItem
from .io.parse.mda import _rejuvenate_universe, _sync_universe
from .io.parse.star import _rehydrate_ensembles
from .props import MolecularNodesObjectProperties
from .ui.node_menu import MN_add_node_menu, draw_node_menus
from . import ui
from .ui.panel import MN_PT_panel, change_style_menu, change_style_node_menu
from .ui import pref

all_classes = (
    ui.CLASSES
    + ops_io
    + [
        MolecularNodesObjectProperties,
        MN_PT_panel,
    ]
    + pref.CLASSES
)

universe_funcs = [_sync_universe, _rejuvenate_universe]


def register():
    # register all of the import operators
    for op in all_classes:
        try:
            bpy.utils.register_class(op)
        except Exception as e:
            print(e)
            pass

    bpy.types.Scene.trajectory_selection_list = bpy.props.CollectionProperty(
        type=TrajectorySelectionItem
    )
    bpy.types.NODE_MT_add.append(MN_add_node_menu)
    bpy.types.Object.mn = bpy.props.PointerProperty(type=MolecularNodesObjectProperties)
    bpy.types.VIEW3D_MT_object_context_menu.prepend(change_style_menu)
    bpy.types.NODE_MT_context_menu.prepend(change_style_node_menu)
    bpy.app.handlers.load_post.append(_rehydrate_ensembles)

    for func in universe_funcs:
        try:
            bpy.app.handlers.load_post.append(func)
        except ValueError as e:
            print(f"Failed to append {func}, error: {e}.")


def unregister():
    # unregister all of the import operator classes
    for op in all_classes:
        try:
            bpy.utils.unregister_class(op)
        except Exception as e:
            print(e)
            pass

    bpy.types.NODE_MT_add.remove(MN_add_node_menu)
    bpy.types.VIEW3D_MT_object_context_menu.remove(change_style_menu)
    bpy.types.NODE_MT_context_menu.remove(change_style_node_menu)
    bpy.app.handlers.load_post.append(_rehydrate_ensembles)

    del bpy.types.Scene.trajectory_selection_list
    try:
        del bpy.types.Object.mn
    except AttributeError:
        print("bpy.types.Object.mn not registered, unable to delete")

    for func in universe_funcs:
        try:
            bpy.app.handlers.load_post.remove(func)
        except ValueError as e:
            print(f"Failed to remove {func}, error: {e}.")


# # # register won't be called when MN is run as a module
bpy.app.handlers.load_post.append(_rejuvenate_universe)
bpy.app.handlers.save_post.append(_sync_universe)
