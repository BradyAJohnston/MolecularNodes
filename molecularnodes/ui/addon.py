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
from bpy.app.handlers import frame_change_pre, load_post, save_post
from bpy.props import CollectionProperty, PointerProperty
from .. import session
from ..handlers import update_entities
from ..utils import add_current_module_to_path
from . import node_menu, ops, panel, pref, props

all_classes = (
    panel.CLASSES
    + ops.CLASSES
    + props.CLASSES
    + pref.CLASSES
    + session.CLASSES
    + node_menu.CLASSES
)


def _test_register():
    try:
        register()
    except Exception as e:
        print(e)
        unregister()
        register()


def register():
    # register all of the import operators
    for op in all_classes:
        try:
            bpy.utils.register_class(op)
        except Exception as e:
            print(e)
            pass
    add_current_module_to_path()
    bpy.types.NODE_MT_add.append(node_menu.add_node_menu)
    bpy.types.VIEW3D_MT_object_context_menu.prepend(panel.pt_object_context)
    bpy.types.NODE_MT_context_menu.prepend(panel.change_style_node_menu)

    save_post.append(session._pickle)
    load_post.append(session._load)
    frame_change_pre.append(update_entities)

    bpy.types.Scene.MNSession = session.MNSession()  # type: ignore
    bpy.types.Object.uuid = props.uuid_property  # type: ignore
    bpy.types.Object.mn = PointerProperty(type=props.MolecularNodesObjectProperties)  # type: ignore
    bpy.types.Scene.mn = PointerProperty(type=props.MolecularNodesSceneProperties)  # type: ignore
    bpy.types.Object.mn_trajectory_selections = CollectionProperty(  # type: ignore
        type=props.TrajectorySelectionItem  # type: ignore
    )


def unregister():
    for op in all_classes:
        try:
            bpy.utils.unregister_class(op)
        except Exception as e:
            print(e)
            pass

    bpy.types.NODE_MT_add.remove(node_menu.add_node_menu)
    bpy.types.VIEW3D_MT_object_context_menu.remove(panel.pt_object_context)
    bpy.types.NODE_MT_context_menu.remove(panel.change_style_node_menu)

    save_post.remove(session._pickle)
    load_post.remove(session._load)
    frame_change_pre.remove(update_entities)
    del bpy.types.Scene.MNSession  # type: ignore
    del bpy.types.Scene.mn  # type: ignore
    del bpy.types.Object.mn  # type: ignore
    del bpy.types.Object.mn_trajectory_selections  # type: ignore
