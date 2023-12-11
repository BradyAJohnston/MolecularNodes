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
    "name"        : "molecularnodes",
    "author"      : "Brady Johnston",
    "description" : "Toolbox for molecular animations in Blender & Geometry Nodes.",
    "blender"     : (4, 0, 0),
    "version"     : (4, 0, 6),
    "location"    : "Scene Properties -> Molecular Nodes",
    "warning"     : "",
    "doc_url"     : "https://bradyajohnston.github.io/MolecularNodes/",
    "tracker_url" : "https://github.com/BradyAJohnston/MolecularNodes/issues",
    "category"    : "Import"
}

import bpy
from . import io
from .io.mda import _rejuvenate_universe, _sync_universe
from .ui.node_menu import MN_add_node_menu
from .io.load import MolecularNodesObjectProperties
from . import auto_load
from .util.utils import template_install

auto_load.init()

universe_funcs = [_sync_universe, _rejuvenate_universe]

def register():
    auto_load.register()
    bpy.types.NODE_MT_add.append(MN_add_node_menu)
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
        auto_load.unregister()
        del bpy.types.Object.mn
        for func in universe_funcs:
            try:
                bpy.app.handlers.load_post.remove(func)
            except ValueError as e:
                print(f"Failed to remove {func}, error: {e}.")
    except RuntimeError as e:
        print("Unable to unregister: {e}")

# if __name__ == "main":
#     register()

# # register won't be called when MN is run as a module
bpy.app.handlers.load_post.append(_rejuvenate_universe)
bpy.app.handlers.save_post.append(_sync_universe)
