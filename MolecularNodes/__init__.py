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
    "description" : "Toolbox for molecular animations in Blender & Geometry Nodes.",
    "blender"     : (3, 5, 0),
    "version"     : (2, 10, 0),
    "location"    : "Scene Properties -> MolecularNodes",
    "warning"     : "",
    "doc_url"     : "https://bradyajohnston.github.io/MolecularNodes/",
    "tracker_url" : "https://github.com/BradyAJohnston/MolecularNodes/issues",
    "category"    : "Import"
}

from . import auto_load
from .mda import _rejuvenate_universe, _sync_universe
from .ui import MN_add_node_menu
import bpy
from . import utils

auto_load.init()

def register():
    auto_load.register()
    bpy.types.NODE_MT_add.append(MN_add_node_menu)
    utils.template_install()

def unregister():
    bpy.types.NODE_MT_add.remove(MN_add_node_menu)
    auto_load.unregister()
    bpy.app.handlers.load_post.remove(_rejuvenate_universe)
    bpy.app.handlers.save_pre.remove(_sync_universe)

# register won't be called when MN is run as a module
bpy.app.handlers.load_post.append(_rejuvenate_universe)
bpy.app.handlers.save_pre.append(_sync_universe)
