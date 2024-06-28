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

from .utils import template_install
from . import auto_load
from .props import MolecularNodesObjectProperties
from .ui.node_menu import MN_add_node_menu, draw_node_menus
from .io.parse.mda import _rejuvenate_universe, _sync_universe
from .io.parse.star import _rehydrate_ensembles
from .ui.panel import change_style_menu, change_style_node_menu
import bpy

bl_info = {
    "name": "molecularnodes",
    "author": "Brady Johnston",
    "description": "Toolbox for molecular animations in Blender & Geometry Nodes.",
    "blender": (4, 1, 0),
    "version": (4, 1, 4),
    "location": "Scene Properties -> Molecular Nodes",
    "warning": "",
    "doc_url": "https://bradyajohnston.github.io/MolecularNodes/",
    "tracker_url": "https://github.com/BradyAJohnston/MolecularNodes/issues",
    "category": "Import",
}

auto_load.init()

universe_funcs = [_sync_universe, _rejuvenate_universe]


def register():
    auto_load.register()
    bpy.types.NODE_MT_add.append(MN_add_node_menu)
    bpy.types.VIEW3D_MT_object_context_menu.prepend(change_style_menu)
    bpy.types.NODE_MT_context_menu.prepend(change_style_node_menu)
    bpy.types.Object.mn = bpy.props.PointerProperty(type=MolecularNodesObjectProperties)
    for func in universe_funcs:
        try:
            bpy.app.handlers.load_post.append(func)
        except ValueError as e:
            print(f"Filaed to append {func}, error: {e}.")
    template_install()


def unregister():
    try:
        bpy.types.NODE_MT_add.remove(MN_add_node_menu)
        bpy.types.VIEW3D_MT_object_context_menu.remove(change_style_menu)
        bpy.types.NODE_MT_context_menu.remove(change_style_node_menu)

        auto_load.unregister()
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
