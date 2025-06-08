"""
Blender operators
"""

import bpy
import MDAnalysis as mda
from bpy.props import (  # type: ignore
    EnumProperty,
    StringProperty,
)
from ..session import get_session
from ..ui.style import STYLE_ITEMS


class MDA_OT_Add_Universe(bpy.types.Operator):
    """
    Operator to add universe
    """

    bl_idname = "mda.add_universe"
    bl_label = "Add universe"
    bl_description = "Add a new MDAnalysis universe"

    topology: StringProperty(  # type: ignore
        name="Topology",
        description="File path for the toplogy file of the universe",
        subtype="FILE_PATH",
        maxlen=0,
    )
    trajectory: StringProperty(  # type: ignore
        name="Trajectory",
        description="File path for the trajectory file of the universe",
        subtype="FILE_PATH",
        maxlen=0,
    )
    name: StringProperty(  # type: ignore
        name="Name",
        description="Name of the universe",
        default="NewUniverse",
        maxlen=0,
    )
    style: EnumProperty(  # type: ignore
        name="Style",
        description="Style for the universe",
        items=STYLE_ITEMS,
        default="spheres",
    )

    def execute(self, context):
        universe = mda.Universe(self.topology, self.trajectory)
        session = get_session()
        session.MDAVis.universes.add(universe, self.style, self.name)
        return {"FINISHED"}

    def invoke(self, context, event):
        sprop = context.scene.mn.mda
        self.name = "u" + str(sprop.next_index)
        return context.window_manager.invoke_props_dialog(self)


class MDA_OT_Delete_Universe(bpy.types.Operator):
    """
    Operator to delete selected universe
    """

    bl_idname = "mda.delete_universe"
    bl_label = "Delete universe"
    bl_description = "Delete MDAnalysis universe"

    @classmethod
    def poll(cls, context):
        return context.scene.mn.mda.active_index != -1

    def execute(self, context):
        sprop = context.scene.mn.mda
        object_name = sprop.universes[sprop.active_index].object.name
        session = get_session()
        session.MDAVis.universes.delete(object_name)
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_confirm(
            self, event, title="Delete universe?"
        )


class MDA_OT_Frame_Selected_Universe(bpy.types.Operator):
    """
    Operator to frame selected universe
    """

    bl_idname = "mda.frame_selected_universe"
    bl_label = "Frame selected universe"
    bl_description = "Zoom to and frame selected universe"

    index: bpy.props.IntProperty()  # type: ignore

    def execute(self, context):
        context.scene.mn.mda.active_index = self.index
        return {"FINISHED"}


class MDA_OT_Frame_All(bpy.types.Operator):
    """
    Operator to frame all universes
    """

    bl_idname = "mda.frame_all"
    bl_label = "Frame All"
    bl_description = "Frame all universes into view"

    def execute(self, context):
        # TODO: Frame all universes
        return {"FINISHED"}


CLASSES = [
    MDA_OT_Add_Universe,
    MDA_OT_Delete_Universe,
    MDA_OT_Frame_Selected_Universe,
    MDA_OT_Frame_All,
]
