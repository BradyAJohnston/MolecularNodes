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
from bpy.app.handlers import (
    frame_change_pre,
    load_post,
    load_pre,
    persistent,
    render_pre,
    save_post,
)
from bpy.props import CollectionProperty, PointerProperty
from .. import session
from ..handlers import render_pre_handler, update_entities
from ..templates import register_templates_menu, unregister_templates_menu
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

_is_registered = False
_mn_session = None
_mn_annotations = None


def _test_register():
    try:
        register()
    except Exception as e:
        print(e)
        unregister()
        register()


@persistent
def _session_restore_draw_handlers():
    if _mn_session is not None:
        _mn_session.add_draw_handlers()
    return None  # execute only once


def register():
    global _is_registered

    if _is_registered:
        return

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
    load_pre.append(session._remove_draw_handlers)
    frame_change_pre.append(update_entities)
    render_pre.append(render_pre_handler)

    if _mn_session is not None:
        bpy.types.Scene.MNSession = _mn_session
    else:
        bpy.types.Scene.MNSession = session.MNSession()
    bpy.types.Object.uuid = props.uuid_property  # type: ignore
    bpy.types.Object.mn = PointerProperty(type=props.MolecularNodesObjectProperties)  # type: ignore
    bpy.types.Scene.mn = PointerProperty(type=props.MolecularNodesSceneProperties)  # type: ignore
    bpy.types.Object.mn_trajectory_selections = CollectionProperty(  # type: ignore
        type=props.TrajectorySelectionItem  # type: ignore
    )
    # bpy.types.Object.mn_annotations is dynamically created and updated based
    # on different annotation types. It has to be a top level property to avoid
    # AttributeError: '_PropertyDeferred' object has no attribute '...'
    if _mn_annotations is not None:
        # if the extension is being re-registered and the modules are already
        # loaded, reuse the saved value
        bpy.types.Object.mn_annotations = _mn_annotations
    register_templates_menu()
    if _mn_session is not None:
        # register a run once timer to restore draw handlers
        bpy.app.timers.register(_session_restore_draw_handlers, first_interval=0.01)

    _is_registered = True


def unregister():
    global _is_registered
    global _mn_session
    global _mn_annotations

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
    load_pre.remove(session._remove_draw_handlers)
    frame_change_pre.remove(update_entities)
    render_pre.remove(render_pre_handler)

    session._remove_draw_handlers(filepath=None)

    _mn_session = bpy.types.Scene.MNSession
    del bpy.types.Scene.MNSession  # type: ignore
    del bpy.types.Scene.mn  # type: ignore
    del bpy.types.Object.mn  # type: ignore
    del bpy.types.Object.mn_trajectory_selections  # type: ignore
    _mn_annotations = bpy.types.Object.mn_annotations
    del bpy.types.Object.mn_annotations  # type: ignore
    unregister_templates_menu()

    _is_registered = False
